
#ifndef FIXED_SIZE_SERIAL_DENSE_MATRIX_H
#define FIXED_SIZE_SERIAL_DENSE_MATRIX_H

#include <Epetra_LAPACK.h>
#include <Epetra_BLAS.h>
#include <Epetra_SerialDenseMatrix.h>

#include <cstring>
#include <vector>
#include <cmath>
#include <ostream>
#include <iostream>

#include "../drt_lib/drt_dserror.H"


namespace LINALG
{

  /// Compile time error check for fixed size serial dense classes
  /*!

    This namespace contains the compile time error checker used by
    FixedSizeSerialDenseMatrix and FixedSizeSerialDenseSolver.

    \author Martin Kuettler
    \date 09/08
   */
  namespace DENSEERROR
  {

    /// Compile time error definitions
    /*!
      A struct that's used for compile time error checking. The
      template argument is the expression to be checked, if it
      evaluates to true nothing happens. If it is false a compile
      error similiar to
      "'Matrix_dimensions_cannot_be_zero' is not a member of 'LINALG::DENSEERROR::Checker<false>'"
      is generated (with gcc).
      Obviously the test expression must be known at compile time.

      This is the (general) definition that is used for expr==true. It
      defines all known errors als empty static inline functions, so
      that the compiler can optimize them out.
     */
    template <bool expr>
    struct Checker
    {
      static inline void Matrix_dimensions_cannot_be_zero() { };
      static inline void Cannot_call_1D_access_function_on_2D_matrix() { };
      static inline void Cannot_compute_determinant_of_nonsquare_matrix() { };
      static inline void Cannot_compute_inverse_of_nonsquare_matrix() { };
      static inline void Transpose_argument_must_be_N_or_T() { };
      static inline void Matrix_size_in_solver_must_be_square() { };
      static inline void Use_FixedSizeSerialDenseSolver_for_matrices_bigger_than_3x3() { };
    };

    /// Compile time error definitions: missing functions raise errors
    /*!
      This is the specialisation for expr==false. It is empty, so that
      the compiler does not find the functions and raises errors.
     */
    template<>
    struct Checker<false>
    {
    };

  }

  /// Internal templated serial dense functions
  /*!
    This namespace contains inlined templated functions for various
    matrix operations. They need the matrix sizes as template
    arguments and take the matrices as Epetra_SerialDenseMatrix, or as
    double*. Normally there is no need for calling them directly,
    because FixedSizeSerialDenseMatrix provides a more convenient
    interface.

    All double* arguments are considered to be pointers to a matrix stored
    in column-major mode, just like the function A() of Epetra- or
    fixed size serial dense matrices return them.

    \author Martin
    \date 09/08
   */
  namespace DENSEFUNCTIONS
  {

    /*
     * Declaration of the functions taking double*
     *
     */

    /// Multiplication: \e out = \e left[^T]*\e right[^T]
    /*!
      Multiply \e left and \e right and store the result in \e out. This
      function takes up to five template parameters: unsigned ints \c i, \c j
      and \c k and optionally two chars \e trans1 and \e trans2. \c i, \c j and \c k
      denote the sizes of the matrices, \e trans1 and \e trans2 say whether
      \e left and \e right should be transposed. 'T' means the respective matrix
      is to be used transposed, 'N' means it is not (the default).

      \note You must specify either none or both of the char template
      arguments.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template<unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(double* out, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(double* out, const double* left, const double* right);

    /// Multiplication: \e out = \e infac * \e left[^T]*\e right[^T]
    /*!
      Multiply \e left and \e right, scale the result by \e infac and store
      it in \e out. This function takes up to five template parameters:
      unsigned ints \c i, \c j and \c k and optionally two chars \e trans1 and
      \e trans2. \c i, \c j and \c k denote the sizes of the matrices, \e trans1 and
      \e trans2 say whether \e left and \e right should be transposed. 'T'
      means the respective matrix is to be used transposed, 'N' means
      it is not (the default).

      \note You must specify either none or both of the char template
      arguments.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right

      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(double* out, const double infac, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(double* out, const double infac, const double* left, const double* right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right
    /*!
      Scale \e out by \e outfac and add \e left*\e right scaled by \e infac. This
      function takes up to five template parameters: unsigned ints \c i, \c j
      and \c k and optionally two chars \e trans1 and \e trans2. \c i, \c j and \c k
      denote the sizes of the matrices, \e trans1 and \e trans2 say whether
      \e left and \e right should be transposed. 'T' means the respective
      matrix is to be used transposed, 'N' means it is not (the
      default).

      \note You must specify either none or both of the char template
      arguments.

      \param outfac
        scalar to multiply with \e out
      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(const double outfac, double* out, const double infac, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double outfac, double* out, const double infac, const double* left, const double* right);

    /// Invert matrix: \e out = inv(\e in)
    /*!
      Invert the matrix \e in and store the result in \e out. To keep a
      common interface there are two template parameters \c i and \c j, but
      they must be the same number. The sizes of \e in and \e out are
      expected to be (\c i)x(\c j), and they must be square.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c j)
      \param in
        pointer to the matrix to be inverted, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline double invert(double* out, const double* in);

    /// Invert matrix: \e mat = inv(\e mat)
    /*!
      Invert the matrix \e mat in place. To keep a common interface there
      are two template parameters \c i and \c j, but they must be the same
      number. The size of \e mat is expected to be (\c i)x(\c j), and it must
      be square.

      \param mat
        pointer to the matrix to be inverted in place, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline double invert(double* mat);

    /// Compute determinant
    /*!
      Computes and returns the determinant of \e mat. To keep a common
      interface there are two template parameters \c i and \c j, but they
      must be the same number. The size of \e mat is expected to be
      (\c i)x(\c j), and it must be square.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return determinant
     */
    template <unsigned int i, unsigned int j>
    inline double determinant(const double* mat);

    /// Copy: \e out = \e in
    /*!
      Copy \e in to \e out. This function takes two template parameters \c i and
      \c j denoting the sizes of the matrices.

      \param out
        pointer to the result matrix, size (\c i)x(\c j)
      \param in
        pointer to the matrix to be copied, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(double* out, const double* in);

    /// Scaled copy: \e out = \e infac * \e in
    /*!
      Scale \e in by \e infac and store the result in \e out. This function takes two template
      parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        pointer to the result matrix, size (\c i)x(\c j)
      \param infac
        scalar to multiply with \e in
      \param in
        pointer to the matrix to read from, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(double* out, const double infac, const double* in);

    /// Addition: \e out = \e outfac * \e out + \e infac * \e in
    /*!
      Scale \e out by \e outfac and add \e infac * \e in to it. This function
      takes two template parameters \c i and \c j denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        pointer to the result matrix, size (\c i)x(\c j)
      \param infac
        scalar to multiply with \e in
      \param in
        pointer to the matrix to be added, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, double* out, const double infac, const double* in);

    /// Addition: \e out = \e left + \e right
    /*!
      Add \e left and \e right and store the result in \e out. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c j)
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(double* out, const double* left,
                          const double* right);

    /// Addition: \e out = \e leftfac * \e left + \e rightfac * \e right
    /*!
      Add \e left and \e right, scaled by \e leftfac and \e rightfac
      respectively. The result is stored in \e out. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c j)
      \param leftfac
        scalar to multiply with \e left
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param rightfac
        scalar to multiply with \e right
      \param right
        pointer to the second factor, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(double* out, const double leftfac, const double* left, const double rightfac, const double* right);

    /// Addition: \e out = \e outfac * \e out + \e leftfac * \e left + \e rightfac * \e right
    /*!
      Scale \e out by \e outfac and add \e left and \e right, scaled by \e leftfac and \e rightfac
      respectively. The result is stored in \e out. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param outfac
        scalar to multiply \e out with
      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c j)
      \param leftfac
        scalar to multiply with \e left
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param rightfac
        scalar to multiply with \e right
      \param right
        pointer to the second factor, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, double* out, const double leftfac, const double* left,
                          const double rightfac, const double* right);

    /// Scale matrix
    /*!
      Scale \e mat by \e fac. This function takes
      two template parameters \c i and \c j denoting the size of \e mat.

      \param fac
        scalar to multiply with \e mat
      \param mat
        pointer to the matrix, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void scaleMatrix(const double fac, double* mat);

    /// Dot product
    /*!
      Return dot product \e left and \e right. This function
      takes two template parameters \c i and \c j denoting the sizes of the matrices.

      \param left
        pointer to the first matrix, size (\c i)x(\c j)
      \param right
        pointer to the second matrix, size (\c i)x(\c j)
      \return dot product
     */
    template <unsigned int i, unsigned int j>
    inline double dot(const double* left, const double* right);

    /// Set matrix to zero
    /*!
      Set matrix \e mat to zero. This function takes two template
      parameters i and j denoting the size of the matrix.

      This is the same as \e putScalar<\c i, \c j>(0.0, \e mat), but it should be faster.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
     */
    template<unsigned int i, unsigned int j>
    inline void clearMatrix(double* mat);

    /// Fill matrix with scalar value
    /*!
      Set every number in \e mat to \e scalar. This function takes two template
      parameters \c i and \c j denoting the size of the matrix.

      \param scalar
        scalar value to be set
      \param mat
        pointer to the matrix, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void putScalar(const double scalar, double* mat);

    /// Calculate absolut values of a matrix
    /*!
      Fill \e out with the absolute values from \e in. This function takes two
      template parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        pointer to the matrix to be set, size (\c i)x(\c j)
      \param in
        pointer to the matrix the values are read from, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void abs(double* out, const double* in);

    /// Calculate reciprocal values of a matrix
    /*!
      Fill \e out with the reciprocal of the values from \e in. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param out
        pointer to the matrix to be set, size (\c i)x(\c j)
      \param in
        pointer to the matrix the values are read from, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void reciprocal(double* out, const double* in);

    /// 1-norm
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return 1-norm of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double norm1(const double* mat);

    /// 2-norm (Euclidean norm)
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return 2-norm of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double norm2(const double* mat);

    /// Inf-norm
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return inf-norm of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double normInf(const double* mat);

    /// Minimum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return minimum value of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double minValue(const double* mat);

    /// Maximum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return maximum value of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double maxValue(const double* mat);

    /// Mean value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return mean value of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double meanValue(const double* mat);

    /*
     * Declaration of the functions taking Epetra_SerialDenseMatrix
     *
     */


    /// Multiplication: \e out = \e left[^T]*\e right[^T]
    /*!
      Multiply \e left and \e right and store the result in \e out. This
      function takes up to five template parameters: unsigned ints \c i, \c j
      and \c k and optionally two chars \e trans1 and \e trans2. \c i, \c j and \c k
      denote the sizes of the matrices, \e trans1 and \e trans2 say whether
      \e left and \e right should be transposed. 'T' means the respective matrix
      is to be used transposed, 'N' means it is not (the default).

      \note You must specify either none or both of the char template
      arguments.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right);

    template<unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right);

    /// Multiplication: \e out = \e infac * \e left[^T]*\e right[^T]
    /*!
      Multiply \e left and \e right, scale the result by \e infac and store
      it in \e out. This function takes up to five template parameters:
      unsigned ints \c i, \c j and \c k and optionally two chars \e trans1 and
      \e trans2. \c i, \c j and \c k denote the sizes of the matrices, \e trans1 and
      \e trans2 say whether \e left and \e right should be transposed. 'T'
      means the respective matrix is to be used transposed, 'N' means
      it is not (the default).

      \note You must specify either none or both of the char template
      arguments.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right);

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right
    /*!
      Scale \e out by \e outfac and add \e left*\e right scaled by \e infac. This
      function takes up to five template parameters: unsigned ints \c i, \c j
      and \c k and optionally two chars \e trans1 and \e trans2. \c i, \c j and \c k
      denote the sizes of the matrices, \e trans1 and \e trans2 say whether
      \e left and \e right should be transposed. 'T' means the respective
      matrix is to be used transposed, 'N' means it is not (the
      default).

      \note You must specify either none or both of the char template
      arguments.

      \param outfac
        scalar to multiply with \e out
      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(const double outfac, Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right);

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double outfac, Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right);

    /// Invert matrix: \e out = inv(\e in)
    /*!
      Invert the matrix \e in and store the result in \e out. To keep a
      common interface there are two template parameters \c i and \c j, but
      they must be the same number. The sizes of \e in and \e out are
      expected to be (\c i)x(\c j), and they must be square.

      \note This function only works for matrices with sizes up to
      3x3. For larger matrices use the FixedSizeSerialDenseSolver.

      \param out
        matrix the result should be stored in, size (\c i)x(\c j)
      \param in
        matrix to be inverted, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline double invert(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in);

    /// Invert matrix: \e mat = inv(\e mat)
    /*!
      Invert the matrix \e mat in place. To keep a common interface there
      are two template parameters \c i and \c j, but they must be the same
      number. The size of \e mat is expected to be (\c i)x(\c j), and it must
      be square.

      \note This function only works for matrices with sizes up to
      3x3. For larger matrices use the FixedSizeSerialDenseSolver.

      \param mat
        matrix to be inverted in place, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline double invert(Epetra_SerialDenseMatrix& mat);

    /// Compute determinant
    /*!
      Computes and returns the determinant of \e mat. To keep a common
      interface there are two template parameters \c i and \c j, but they
      must be the same number. The size of \e mat is expected to be
      (\c i)x(\c j), and it must be square.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return determinant
     */

    template <unsigned int i, unsigned int j>
    inline double determinant(const Epetra_SerialDenseMatrix& mat);

    /// Copy: \e out = \e in
    /*!
      Copy \e in to \e out. This function takes two template parameters \c i and
      \c j denoting the sizes of the matrices.

      \param out
        result matrix, size (\c i)x(\c j)
      \param in
        matrix to be copied, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in);

    /// Scaled copy: \e out = \e infac * \e in
    /*!
      Scale \e in by \e infac and store the result in \e out. This function takes two template
      parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        result matrix, size (\c i)x(\c j)
      \param infac
        scalar to multiply with \e in
      \param in
        matrix to read from, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const double infac, const Epetra_SerialDenseMatrix& in);

    /// Addition: \e out = \e outfac * \e out + \e infac * \e in
    /*!
      Scale \e out by \e outfac and add \e infac * \e in to it. This function
      takes two template parameters \c i and \c j denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        result matrix, size (\c i)x(\c j)
      \param infac
        scalar to multiply with \e in
      \param in
        matrix to be added, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, Epetra_SerialDenseMatrix& out, const double infac, const Epetra_SerialDenseMatrix& in);

    /// Addition: \e out = \e left + \e right
    /*!
      Add \e left and \e right and store the result in \e out. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c j)
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left,
                          const Epetra_SerialDenseMatrix& right);

    /// Addition: \e out = \e leftfac * \e left + \e rightfac * \e right
    /*!
      Add \e left and \e right, scaled by \e leftfac and \e rightfac
      respectively. The result is stored in \e out. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c j)
      \param leftfac
        scalar to multiply with \e left
      \param left
        first factor, size (\c i)x(\c j)
      \param rightfac
        scalar to multiply with \e right
      \param right
        second factor, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const double leftfac, const Epetra_SerialDenseMatrix& left,
                          const double rightfac, const Epetra_SerialDenseMatrix& right);

    /// Addition: \e out = \e outfac * \e out + \e leftfac * \e left + \e rightfac * \e right
    /*!
      Scale \e out by \e outfac and add \e left and \e right, scaled by \e leftfac and \e rightfac
      respectively. The result is stored in \e out. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param outfac
        scalar to multiply \e out with
      \param out
        matrix the result should be stored in, size (\c i)x(\c j)
      \param leftfac
        scalar to multiply with \e left
      \param left
        first factor, size (\c i)x(\c j)
      \param rightfac
        scalar to multiply with \e right
      \param right
        second factor, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, Epetra_SerialDenseMatrix& out, const double leftfac,
                          const Epetra_SerialDenseMatrix& left, const double rightfac, const Epetra_SerialDenseMatrix& right);

    /// Scale matrix
    /*!
      Scale \e mat by \e scalar. This function takes
      two template parameters \c i and \c j denoting the size of \e mat.

      \param scalar
        scalar to multiply with \e mat
      \param mat
        matrix, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void scaleMatrix(const double scalar, Epetra_SerialDenseMatrix& mat)
    { scaleMatrix<i,j>(scalar, mat.A()); }

    /// Dot product
    /*!
      Return dot product \e left and \e right. This function
      takes two template parameters \c i and \c j denoting the sizes of the matrices.

      \param left
        first matrix, size (\c i)x(\c j)
      \param right
        second matrix, size (\c i)x(\c j)
      \return dot product
     */
    template <unsigned int i, unsigned int j>
    inline double dot(const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    { return dot<i,j>(left.A(), right.A()); }

    /// Set matrix to zero
    /*!
      Set matrix \e mat to zero. This function takes two template
      parameters i and j denoting the size of the matrix.

      This is the same as \e putScalar<\c i,\c j>(0.0, \e mat), but it should be faster.

      \param mat
        matrix, size (\c i)x(\c j)
     */
    template<unsigned int i, unsigned int j>
    inline void clearMatrix(Epetra_SerialDenseMatrix& mat)
    { clearMatrix<i,j>(mat.A()); }

    /// Fill matrix with scalar value
    /*!
      Set every number in \e mat to \e scalar. This function takes two template
      parameters \c i and \c j denoting the size of the matrix.

      \param scalar
        scalar value to be set
      \param mat
        matrix, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void putScalar(const double scalar, Epetra_SerialDenseMatrix& mat)
    { putScalar<i,j>(scalar, mat.A()); }

    /// Calculate absolut values of a matrix
    /*!
      Fill \e out with the absolute values from \e in. This function takes two
      template parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        matrix to be set, size (\c i)x(\c j)
      \param in
        matrix the values are read from, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void abs(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in);

    /// Calculate reciprocal values of a matrix
    /*!
      Fill \e out with the reciprocal of the values from \e in. This
      function takes two template parameters \c i and \c j denoting the
      sizes of the matrices.

      \param out
        matrix to be set, size (\c i)x(\c j)
      \param in
        matrix the values are read from, size (\c i)x(\c j)
     */
    template <unsigned int i, unsigned int j>
    inline void reciprocal(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in);

    /// 1-norm
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return 1-norm of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double norm1(const Epetra_SerialDenseMatrix& mat)
    { return norm1<i,j>(mat.A()); }

    /// 2-norm (Euclidean norm)
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return 2-norm of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double norm2(const Epetra_SerialDenseMatrix& mat)
    { return norm2<i,j>(mat.A()); }

    /// Inf-norm
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return inf-norm of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double normInf(const Epetra_SerialDenseMatrix& mat)
    { return normInf<i,j>(mat.A()); }

    /// Minimum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return minimum value of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double minValue(const Epetra_SerialDenseMatrix& mat)
    { return minValue<i,j>(mat.A()); }

    /// Maximum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return maximum value of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double maxValue(const Epetra_SerialDenseMatrix& mat)
    { return maxValue<i,j>(mat.A()); }

    /// Mean value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return mean value of \e mat
     */
    template <unsigned int i,  unsigned int j>
    inline double meanValue(const Epetra_SerialDenseMatrix& mat)
    { return meanValue<i,j>(mat.A()); }

    /*
     * Definitions of the functions taking double*
     *
     */

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(double* out, const double* left, const double* right)
    {
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          *out = (*left)*(*right);
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            left += i;
            ++right;
            *out += (*left)*(*right);
          }
          ++out;
          left -= i*(j-1) - 1;
          right -= j-1;
        }
        left -= i;
        right += j;
      }
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline
    void multiply(double* out, const double* left, const double* right)
    {
      DENSEERROR::Checker<(trans1 == 'N' or trans1 == 'T') and (trans2 == 'N' or trans2 == 'T')>::Transpose_argument_must_be_N_or_T();
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          *out = (*left)*(*right);
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            switch(trans1) {
            case 'N': left += i; break;
            case 'T':    ++left; break;
            }
            switch(trans2) {
            case 'N':    ++right; break;
            case 'T': right += k; break;
            }
            *out += (*left)*(*right);
          }
          ++out;
          switch(trans1) {
          case 'N': left -= i*(j-1) - 1; break;
          case 'T':              ++left; break;
          }
          switch(trans2) {
          case 'N':     right -= j-1; break;
          case 'T': right -= (j-1)*k; break;
          }
        }
        switch(trans1) {
        case 'N':   left -= i; break;
        case 'T': left -= i*j; break;
        }
        switch(trans2) {
        case 'N': right += j; break;
        case 'T':    ++right; break;
        }
      }
    }

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(double* out, double infac, const double* left, const double* right)
    {
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          *out = (*left)*(*right);
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            left += i;
            ++right;
            *out += (*left)*(*right);
          }
          *out *= infac;
          ++out;
          left -= i*(j-1) - 1;
          right -= j-1;
        }
        left -= i;
        right += j;
      }
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline
    void multiply(double* out, const double infac, const double* left, const double* right)
    {
      DENSEERROR::Checker<(trans1 == 'N' or trans1 == 'T') and (trans2 == 'N' or trans2 == 'T')>::Transpose_argument_must_be_N_or_T();
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          *out = (*left)*(*right);
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            switch(trans1) {
            case 'N': left += i; break;
            case 'T':    ++left; break;
            }
            switch(trans2) {
            case 'N':    ++right; break;
            case 'T': right += k; break;
            }
            *out += (*left)*(*right);
          }
          *out *= infac;
          ++out;
          switch(trans1) {
          case 'N': left -= i*(j-1) - 1; break;
          case 'T':              ++left; break;
          }
          switch(trans2) {
          case 'N':     right -= j-1; break;
          case 'T': right -= (j-1)*k; break;
          }
        }
        switch(trans1) {
        case 'N':   left -= i; break;
        case 'T': left -= i*j; break;
        }
        switch(trans2) {
        case 'N': right += j; break;
        case 'T':    ++right; break;
        }
      }
    }

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double outfac, double* out, double infac, const double* left, const double* right)
    {
      double tmp;
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          tmp = (*left)*(*right);
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            left += i;
            ++right;
            tmp += (*left)*(*right);
          }
          *out *= outfac;
          *out += tmp * infac;
          ++out;
          left -= i*(j-1) - 1;
          right -= j-1;
        }
        left -= i;
        right += j;
      }
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline
    void multiply(const double outfac, double* out, const double infac, const double* left, const double* right)
    {
      DENSEERROR::Checker<(trans1 == 'N' or trans1 == 'T') and (trans2 == 'N' or trans2 == 'T')>::Transpose_argument_must_be_N_or_T();
      double tmp;
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          tmp = (*left)*(*right);
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            switch(trans1) {
            case 'N': left += i; break;
            case 'T':    ++left; break;
            }
            switch(trans2) {
            case 'N':    ++right; break;
            case 'T': right += k; break;
            }
            tmp += (*left)*(*right);
          }
          *out *= outfac;
          *out += tmp * infac;
          ++out;
          switch(trans1) {
          case 'N': left -= i*(j-1) - 1; break;
          case 'T':              ++left; break;
          }
          switch(trans2) {
          case 'N':     right -= j-1; break;
          case 'T': right -= (j-1)*k; break;
          }
        }
        switch(trans1) {
        case 'N':   left -= i; break;
        case 'T': left -= i*j; break;
        }
        switch(trans2) {
        case 'N': right += j; break;
        case 'T':    ++right; break;
        }
      }
    }

    template <unsigned int i, unsigned int j>
    inline double invert(double* out, const double* in)
    {
      DENSEERROR::Checker<i==j>::Cannot_compute_inverse_of_nonsquare_matrix();
      DENSEERROR::Checker<(i < 4)>::Use_FixedSizeSerialDenseSolver_for_matrices_bigger_than_3x3();
      return 0.0;
    }

    template<>
    inline double invert<1,1>(double* out, const double* in)
    {
      const double det = in[0];
      if (det == 0.0) dserror("Determinant of 1x1 matrix is zero");
      out[0] = 1.0/in[0];
      return det;
    }

    template<>
    inline double invert<2,2>(double* out, const double* in)
    {
      const double det = in[0]*in[1+1*2] - in[1]*in[1*2];
      if (det == 0.0) dserror("Determinant of 2x2 matrix is zero");
      const double invdet = 1.0/det;
      out[0    ] =  invdet*in[1+1*2];
      out[1    ] = -invdet*in[1    ];
      out[  1*2] = -invdet*in[  1*2];
      out[1+1*2] =  invdet*in[0    ];
      return det;
    }


    template<>
    inline double invert<3,3>(double* out, const double* in)
    {
      out[0] = in[1+1*3]*in[2+2*3] - in[2+1*3]*in[1+2*3];
      out[1] = in[2    ]*in[1+2*3] - in[1    ]*in[2+2*3];
      out[2] = in[1    ]*in[2+1*3] - in[2    ]*in[1+1*3];
      const double det = in[0] * out[0] + in[1*3] * out[1] + in[2*3] * out[2];
      //const double det = in[0]*in[1+3*1]*in[2+3*2] +
      //                   in[0+3*1]*in[1+3*2]*in[2+3*0] +
      //                   in[0+3*2]*in[1+3*0]*in[2+3*1] -
      //                   in[0+3*2]*in[1+3*1]*in[2+3*0] -
      //                   in[0+3*0]*in[1+3*2]*in[2+3*1] -
      //                   in[0+3*1]*in[1+3*0]*in[2+3*2];
      if (det == 0.0) dserror("Determinant of 3x3 matrix is zero");
      const double invdet = 1.0/det;
      out[0    ] *= invdet;
      out[1    ] *= invdet;
      out[2    ] *= invdet;
      out[  1*3] = invdet * (in[2+1*3]*in[  2*3] - in[  1*3]*in[2+2*3]);
      out[1+1*3] = invdet * (in[0    ]*in[2+2*3] - in[2    ]*in[  2*3]);
      out[2+1*3] = invdet * (in[2    ]*in[  1*3] - in[0    ]*in[2+1*3]);
      out[  2*3] = invdet * (in[  1*3]*in[1+2*3] - in[1+1*3]*in[  2*3]);
      out[1+2*3] = invdet * (in[1    ]*in[  2*3] - in[0    ]*in[1+2*3]);
      out[2+2*3] = invdet * (in[0    ]*in[1+1*3] - in[1    ]*in[  1*3]);
      return det;
    }

    template <unsigned int i, unsigned int j>
    inline double invert(double* mat)
    {
      DENSEERROR::Checker<i==j>::Cannot_compute_inverse_of_nonsquare_matrix();
      DENSEERROR::Checker<(i < 4)>::Use_FixedSizeSerialDenseSolver_for_matrices_bigger_than_3x3();
      return 0.0;
    }

    template <>
    inline double invert<1,1>(double* mat)
    {
      const double det = mat[0];
      if (det == 0.0) dserror("Determinant of 1x1 matrix is zero");
      mat[0] = 1.0/mat[0];
      return det;
    }

    template<>
    inline double invert<2,2>(double* mat)
    {
      double tmp;
      const double det = mat[0]*mat[1+1*2] - mat[1]*mat[1*2];
      if (det == 0.0) dserror("Determinant of 2x2 matrix is zero");
      const double invdet = 1.0/det;
      tmp = mat[0    ]; mat[0    ] = invdet*mat[1+1*2]; mat[1+1*2] = invdet*tmp;
      mat[1    ] *= -invdet;
      mat[  1*2] *= -invdet;
      return det;
    }


    template<>
    inline double invert<3,3>(double* mat)
    {
      const double tmp00 = mat[1+1*3]*mat[2+2*3] - mat[2+1*3]*mat[1+2*3];
      const double tmp10 = mat[2    ]*mat[1+2*3] - mat[1    ]*mat[2+2*3];
      const double tmp20 = mat[1    ]*mat[2+1*3] - mat[2    ]*mat[1+1*3];
      const double det = mat[0] * tmp00 + mat[1*3] * tmp10 + mat[2*3] * tmp20;
      //const double det = mat[0+3*0]*mat[1+3*1]*mat[2+3*2] +
      //                    mat[0+3*1]*mat[1+3*2]*mat[2+3*0] +
      //                    mat[0+3*2]*mat[1+3*0]*mat[2+3*1] -
      //                    mat[0+3*2]*mat[1+3*1]*mat[2+3*0] -
      //                    mat[0+3*0]*mat[1+3*2]*mat[2+3*1] -
      //                    mat[0+3*1]*mat[1+3*0]*mat[2+3*2];
      if (det == 0.0) dserror("Determinant of 3x3 matrix is zero");
      const double invdet = 1.0/det;
      const double tmp01 = mat[  1*3];
      const double tmp11 = mat[1+1*3];
      const double tmp12 = mat[1+2*3];
      mat[  1*3] = invdet * (mat[2+1*3] * mat[  2*3] -    tmp01   * mat[2+2*3]);
      mat[1+1*3] = invdet * (mat[0    ] * mat[2+2*3] - mat[2    ] * mat[  2*3]);
      mat[1+2*3] = invdet * (mat[1    ] * mat[  2*3] - mat[0    ] *    tmp12  );
      mat[2+1*3] = invdet * (mat[2    ] *    tmp01   - mat[0    ] * mat[2+1*3]);
      mat[  2*3] = invdet * (   tmp01   *    tmp12   -    tmp11   * mat[  2*3]);
      mat[2+2*3] = invdet * (mat[0    ] *    tmp11   - mat[1    ] *    tmp01  );
      mat[0    ] = invdet * tmp00;
      mat[1    ] = invdet * tmp10;
      mat[2    ] = invdet * tmp20;
      return det;
    }


    template <unsigned int i, unsigned int j>
    inline double determinant(const double* mat)
    {
      DENSEERROR::Checker<i == j>::Cannot_compute_determinant_of_nonsquare_matrix();

      // taken from src/drt_lib/linalg_utils.cpp: LINALG::DeterminantLU,
      // only with minor changes.
      double* tmp = new double[i*j];
      update<i,j>(tmp,mat);
      Epetra_LAPACK lapack;
      std::vector<int> ipiv(j);
      int info;
      lapack.GETRF(i,j,tmp,i,&ipiv[0],&info);
      if (info<0) dserror("Lapack's dgetrf returned %d",info);
      else if (info>0) return 0.0;
      double d = *tmp;
      for (unsigned int c=1; c<j; ++c) d *= tmp[c+i*c];
      // swapping rows of A changes the sign of the determinant, so we have to
      // undo lapack's permutation w.r.t. the determinant
      // note the fortran indexing convention in ipiv
      for (unsigned int c=0; c<j; ++c)
        if (static_cast<unsigned>(ipiv[c])!=c+1) d *= -1.0;
      return d;
    }

    template <>
    inline double determinant<1,1>(const double* mat)
    {
      return *mat;
    }

    template <>
    inline double determinant<2,2>(const double* mat)
    {
      return mat[0    ]*mat[1+1*2] - mat[1    ]*mat[  1*2];
    }

    template <>
    inline double determinant<3,3>(const double* mat)
    {
      return mat[0    ] * (mat[1+1*3]*mat[2+2*3] - mat[2+1*3]*mat[1+2*3]) +
        mat[  1*3] * (mat[2    ]*mat[1+2*3] - mat[1    ]*mat[2+2*3]) +
        mat[  2*3] * (mat[1    ]*mat[2+1*3] - mat[2    ]*mat[1+1*3]);
      //return mat[0]*mat[1+3*1]*mat[2+3*2] +
      //  mat[0+3*1]*mat[1+3*2]*mat[2+3*0] +
      //  mat[0+3*2]*mat[1+3*0]*mat[2+3*1] -
      //  mat[0+3*2]*mat[1+3*1]*mat[2+3*0] -
      //  mat[0+3*0]*mat[1+3*2]*mat[2+3*1] -
      //  mat[0+3*1]*mat[1+3*0]*mat[2+3*2];
    }

    /* add matrices */

    template <unsigned int i, unsigned int j>
    inline
    void update(double* out, const double* left, const double* right)
    {
      *out = *left + *right;
      for (unsigned int c = 1; c < i*j; ++c)
        *(++out) = *(++left) + *(++right);
    }

    template <unsigned int i, unsigned int j>
    inline
    void update(double* out, const double leftfac, const double* left, const double rightfac, const double* right)
    {
      *out = leftfac*(*left) + rightfac*(*right);
      for (unsigned int c = 1; c < i*j; ++c)
        *(++out) = leftfac*(*(++left)) + rightfac*(*(++right));
    }

    template <unsigned int i, unsigned int j>
    inline
    void update(const double outfac, double* out, const double leftfac, const double* left,
                   const double rightfac, const double* right)
    {
      *out *= outfac;
      *out += leftfac*(*left) + rightfac*(*right);
      for (unsigned int c = 1; c < i*j; ++c) {
        *(++out) *= outfac;
        *out += leftfac*(*(++left)) + rightfac*(*(++right));
      }
    }

    template <unsigned int i, unsigned int j>
    inline void update(double* out, const double* in)
    {
      std::memcpy(out, in, i*j*sizeof(double));
    }

    template <unsigned int i, unsigned int j>
    inline void update(double* out, const double infac, const double* in)
    {
      *out = infac*(*in);
      for (unsigned int c = 1; c < i*j; ++c)
        *(++out) = infac*(*(++in));
    }

    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, double* out, const double infac, const double* in)
    {
      *out *= outfac;
      *out += infac*(*in);
      for (unsigned int c = 1; c < i*j; ++c){
        *(++out) *= outfac;
        *out += infac*(*(++in));
      }
    }

    template <unsigned int i, unsigned int j>
    inline
    void scaleMatrix(const double scalar, double* mat)
    {
      *mat *= scalar;
      for (unsigned int c = 1; c < i*j; ++c)
        *(++mat) *= scalar;
    }

    template <unsigned int i, unsigned int j>
    inline double dot(const double* left, const double* right)
    {
      double res = (*left)*(*right);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++left; ++right;
        res += (*left)*(*right);
      }
      return res;
    }

    template <unsigned int i, unsigned int j>
    inline void clearMatrix(double* mat)
    {
      std::memset(mat,0,i*j*sizeof(double));
    }


    template <unsigned int i, unsigned int j>
    inline void putScalar(const double scalar, double* mat)
    {
      *mat = scalar;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        *mat = scalar;
      }
    }

    template <unsigned int i, unsigned int j>
    inline void abs(double* out, const double* in)
    {
      *out = *in >= 0 ? *in : -*in;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++out; ++in;
        *out = *in >= 0 ? *in : -*in;
      }
    }

    template <unsigned int i, unsigned int j>
    inline void reciprocal(double* out, const double* in)
    {
      *out = 1.0/(*in);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++out; ++in;
        *out = 1.0/(*in);
      }
    }

    template <unsigned int i,  unsigned int j>
    inline double norm1(const double* mat)
    {
      double result = *mat >= 0 ? *mat : -(*mat);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        result += *mat >= 0 ? *mat : -(*mat);
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double norm2(const double* mat)
    {
      double result = (*mat)*(*mat);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        result += (*mat)*(*mat);
      }
      return std::sqrt(result);
    }

    template <unsigned int i,  unsigned int j>
    inline double normInf(const double* mat)
    {
      double result = *mat > 0 ? *mat : -*mat;
      double tmp;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        tmp = *mat >= 0 ? *mat : -*mat;
        if (tmp > result) result = tmp;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double minValue(const double* mat)
    {
      double result = *mat;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        if (*mat < result) result = *mat;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double maxValue(const double* mat)
    {
      double result = *mat;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        if (*mat > result) result = *mat;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double meanValue(const double* mat)
    {
      double result = *mat;
      for (unsigned int c = 1; c < i*j; ++c)
        result += *(++mat);
      return result/(i*j);
    }

#if 0
    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(double* out, const double* left, const double* right)
    {
      *out = (*right)/(*left);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++out; ++right; ++left;
        *out = (*right)/(*left);
      }
    }

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(const double outfac, double* out, const double prodfac, const double* left, const double* right)
    {
      *out = *out*outfac + prodfac*(*right)/(*left);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++out; ++right; ++left;
        *out = *out*outfac + prodfac*(*right)/(*left);
      }
    }
#endif

     /*
     * Definitions of the functions taking Epetra_SerialDenseMatrix
     *
     */

    template<unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
        dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      multiply<i,j,k>(out.A(), left.A(), right.A());
    }


    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      switch (trans1) {
      case 'N':
        switch (trans2) {
        case 'N':
          if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.M() or out.N() != right.M() or left.N() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      case 'T':
        switch (trans2) {
        case 'N':
          if (out.M() != left.N() or out.N() != right.N() or left.M() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.N() or out.N() != right.M() or left.M() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      }
#endif
      multiply<i,j,k,trans1,trans2>(out.A(), left.A(), right.A());
    }

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
        dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      multiply<i,j,k>(out.A(), infac, left.A(), right.A());
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      switch (trans1) {
      case 'N':
        switch (trans2) {
        case 'N':
          if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.M() or out.N() != right.M() or left.N() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      case 'T':
        switch (trans2) {
        case 'N':
          if (out.M() != left.N() or out.N() != right.N() or left.M() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.N() or out.N() != right.M() or left.M() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      }
#endif
      multiply<i,j,k,trans1,trans2>(out.A(), infac, left.A(), right.A());
    }


    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double outfac, Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
        dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      multiply<i,j,k>(outfac, out.A(), infac, left.A(), right.A());
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(const double outfac, Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      switch (trans1) {
      case 'N':
        switch (trans2) {
        case 'N':
          if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.M() or out.N() != right.M() or left.N() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      case 'T':
        switch (trans2) {
        case 'N':
          if (out.M() != left.N() or out.N() != right.N() or left.M() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.N() or out.N() != right.M() or left.M() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      }
#endif
      multiply<i,j,k,trans1,trans2>(outfac, out.A(), infac, left.A(), right.A());
    }

    template <unsigned int i, unsigned int j>
    inline double invert(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != out.N() or in.M() != in.N() or out.M() != in.M())
        dserror("Invalid matrix sizes for inversion, (%i,%i) = inv( (%i,%i) )",out.M(),out.N(),in.M(),in.N());
#endif
      return invert<i,j>(out.A(), in.A());
    }

    template <unsigned int i, unsigned int j>
    inline double invert(Epetra_SerialDenseMatrix& mat)
    {
#ifdef DEBUG
      if (mat.M() != mat.N())
        dserror("Invalid matrix sizes for inversion, inv( (%i,%i) )",mat.M(),mat.N());
#endif
      return invert<i,j>(mat.A());
    }

    template <unsigned int i, unsigned int j>
    inline double determinant(const Epetra_SerialDenseMatrix& mat)
    {
#ifdef DEBUG
      if (mat.M() != mat.N())
        dserror("Invalid matrix sizes for determinant, inv( (%i,%i) )",mat.M(),mat.N());
#endif
      return determinant<i,j>(mat.A());
    }

    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left,
                       const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or left.M() != right.M() or out.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      update<i,j>(out.A(), left.A(), right.A());
    }

    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const double leftfac, const Epetra_SerialDenseMatrix& left,
                       const double rightfac, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or left.M() != right.M() or out.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      update<i,j>(out.A(), leftfac, left.A(), rightfac, right.A());
    }

    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, Epetra_SerialDenseMatrix& out, const double leftfac,
                       const Epetra_SerialDenseMatrix& left, const double rightfac, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or left.M() != right.M() or out.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      update<i,j>(outfac, out.A(), leftfac, left.A(), rightfac, right.A());
    }

    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != in.M() or out.N() != in.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)",
                out.M(),out.N(), in.M(), in.N());
#endif
      update<i,j>(out.A(), in.A());
    }

    template <unsigned int i, unsigned int j>
    inline void update(Epetra_SerialDenseMatrix& out, const double infac, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != in.M() or out.N() != in.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)",
                out.M(),out.N(), in.M(), in.N());
#endif
      update<i,j>(out.A(), infac, in.A());
    }


    template <unsigned int i, unsigned int j>
    inline void update(const double outfac, Epetra_SerialDenseMatrix& out, const double infac, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != in.M() or out.N() != in.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)",
                out.M(),out.N(), in.M(), in.N());
#endif
      update<i,j>(outfac, out.A(), infac, in.A());
    }

    template <unsigned int i, unsigned int j>
    inline void abs(Epetra_SerialDenseMatrix& dest, const Epetra_SerialDenseMatrix& src)
    {
#ifdef DEBUG
      if (dest.M() != dest.M() or src.N() != src.N())
        dserror("Invalid matrix sizes for abs, (%i,%i) = abs( (%i,%i) )",
                dest.M(),dest.N(), src.M(), src.N());
#endif
      abs<i,j>(dest.A(), src.A());
    }

    template <unsigned int i, unsigned int j>
    inline void reciprocal(Epetra_SerialDenseMatrix& dest, const Epetra_SerialDenseMatrix& src)
    {
#ifdef DEBUG
      if (dest.M() != dest.M() or src.N() != src.N())
        dserror("Invalid matrix sizes for reciprocal, (%i,%i) = reciprocal( (%i,%i) )",
                dest.M(),dest.N(), src.M(), src.N());
#endif
      reciprocal<i,j>(dest.A(), src.A());
    }

#if 0
    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(Epetra_SerialDenseMatrix& dest, const Epetra_SerialDenseMatrix& left,
                                   const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (dest.M() != left.M() or left.M() != right.M() or dest.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for reciprocal multiplication, (%i,%i) = reciprocalMultiply( (%i,%i), (%i,%i) )",
                dest.M(),dest.N(), left.M(), left.N(), right.M(), right.N());
#endif
      reciprocalMultiply<i,j>(dest.A(), left.A(), right.A());
    }

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(const double destfac, Epetra_SerialDenseMatrix& dest, const double prodfac,
                                   const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (dest.M() != left.M() or left.M() != right.M() or dest.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for reciprocal multiplication, (%i,%i) = reciprocalMultiply( (%i,%i), (%i,%i) )",
                dest.M(),dest.N(), left.M(), left.N(), right.M(), right.N());
#endif
      reciprocalMultiply<i,j>(destfac, dest.A(), prodfac, left.A(), right.A());
    }
#endif
  }


/// Serial dense matrix with templated dimensions for maximum speed
/*!
  A serial dense matrix with templated dimensions that is supposed to
  be fast and lightweight.
  The double-array is allocated on the heap and stored in
  column-major order, just like in Epetra_SerialDenseMatrix.

  The interface is based on that of Epetra_SerialDenseMatrix and
  Epetra_MultiVector. The whole View/Copy thing works a little
  different, though. See the appropriate functions for details.

  There is no operator[]. It behaves differently in
  Epetra_SerialDenseMatrix and Epetra_SerialDenseVector, and is not
  needed in either of them.

  \author Martin
  \date 09/08
 */
template <unsigned int rows, unsigned int cols /* =1  ?*/ >
class FixedSizeSerialDenseMatrix
{
private:
  /// the pointer holding the data
  double* data_;
  /// whether we are a view to some other matrix
  bool isview_;

public:
  /// Default constructor
  /*!
    Constructs a new FixedSizeSerialDenseMatrix and allocates the
    memory. If \e setzero==true it is filled with zeros, otherwise it
    is left uninitialized.

    \param setzero
      whether matrix should be initialised to zero
   */
  FixedSizeSerialDenseMatrix(bool setzero = false);

  /// Constructor
  /*!
    Constructs a new FixedSizeSerialDenseMatrix from data \e d. If
    \e view==false (the default) the data is copied, otherwise a view to
    it is constructed.

    \param d
      pointer to data
    \param view
      whether the data is to be viewed or copied
   */
  explicit FixedSizeSerialDenseMatrix(double* d, bool view=false);

  /// Constructor
  /*!
    Constructs a new FixedSizeSerialDenseMatrix from \e source. If
    \e view==false the data is copied, otherwise a view to
    it is constructed.

    When both an Epetra and a fixed size version of a matrix is needed
    I recommend constructing an Epetra matrix and having a fixed size
    view onto it. That's because Epetra-Views behave differently than
    normal Epetra matrices in some ways, which can lead to tricky bugs.

    \param source
      matrix to take data from
    \param view
      whether the data is to be viewed or copied
   */
  FixedSizeSerialDenseMatrix(FixedSizeSerialDenseMatrix<rows,cols>& source, bool view);

  /// Copy constructor
  /*!
    Constructs a new FixedSizeSerialDenseMatrix from source. Unlike
    the Epetra_SerialDenseMatrix copy constructor this one *always*
    copies the data, even when source is a view.

    \param source
      matrix to copy
   */
  FixedSizeSerialDenseMatrix(const FixedSizeSerialDenseMatrix<rows,cols>& source);

  /// Deconstructor
  ~FixedSizeSerialDenseMatrix();

  /// Return the double* holding the data.
  inline const double* A() const { return data_; }
  /// Return the double* holding the data.
  inline double* A() { return data_; }
  /// Return the number of rows
  inline unsigned int Rows() const { return rows; }
  /// Return the number of columns
  inline unsigned int Columns() const { return cols; }
  /// Return the number of columns
  inline unsigned int Cols() const { return cols; }
  /// Return the number of rows
  inline unsigned int M() const { return rows; }
  /// Return the number of columns
  inline unsigned int N() const { return cols; }
  /// Check whether the matrix is initialized
  /*!
    You cannot test whether the matrix is empty using M() and N(),
    for they will always return the templated size. Instead this
    function can be used, it tests whether the data pointer is not
    NULL.

    \note To actually get a matrix for which IsInitialized() returns
    false you must construct a view to NULL, because the default
    constructor already allocates memory.
   */
  inline bool IsInitialized() const { return data_!=NULL; }

  /// Set view
  /*!
    Set this matrix to be a view to \e data.

    \param data
      memory to be viewed
   */
  void SetView(double* data);

  /// Set view
  /*!
    Set this matrix to be a view to \e source.

    \param source
      matrix to be viewed
   */
  void SetView(FixedSizeSerialDenseMatrix<rows, cols>& source);

  /// Set copy
  /*!
    Set this matrix to be a copy of \e data. The difference to Update(\e data)
    is that this funcion will allocate it's own memory when it was a
    view before, Update would copy the data into the view.

    \param data
      memory to copy
   */
  void SetCopy(double* data);

  /// Set copy
  /*!
    Set this matrix to be a copy of source. Only the double array
    will be copied, the \e isview_ flag is ignored (this is equivalent to SetCopy(source.A()).
    The difference to Update(\e source) is that this funcion will allocate
    it's own memory when it was a view before, Update would copy the
    data into the view.

    \param source
      matrix to copy from
   */
  void SetCopy(FixedSizeSerialDenseMatrix<rows, cols>& source);


  /// Calculate determinant
  /*!
    \return determinant
   */
  inline double Determinant() const;

  /// Invert in place
  /*!
    Invert this matrix in place.

    \return determinant of matrix before inversion
   */
  inline double Invert();

  /// Invert matrix
  /*!
    Invert matrix \e other and store the result in \e this.

    \param other
      matrix to be inverted
    \return determinant of \e other
   */
  inline double Invert(const FixedSizeSerialDenseMatrix<rows, cols>& other);

  /// Set to zero
  /*!
    Sets every value in this matrix to zero. This is equivalent to
    PutScalar(0.0), but it should be faster.
   */
  inline void Clear();

  // Epetra-style Functions
  /// Fill with scalar
  /*!
    Sets every value in this matrix to \e scalar.

    \param scalar
      value to fill matrix with
   */
  inline void PutScalar(const double scalar);

  /// Dot product
  /*!
    Return the dot product of \e this and \e other.

    \param other
      second factor
    \return dot product
   */
  inline double Dot(const FixedSizeSerialDenseMatrix<rows, cols>& other);

  /// Compute absolute value
  /*!
    Fill this matrix with the absolute value of the numbers in \e other.

    \param other
      matrix to read values from
   */
  inline void Abs(const FixedSizeSerialDenseMatrix<rows, cols>& other);

  /// Compute reciprocal value
  /*!
    Fill this matrix with the reciprocal value of the numbers in \e other.

    \param other
      matrix to read values from
   */
  inline void Reciprocal(const FixedSizeSerialDenseMatrix<rows, cols>& other);

  /// Scale
  /*!
    Scale matrix with \e scalar.

    \param scalar
      scaling factor
   */
  inline void Scale(const double scalar);

  /// Copy: \e this = \e other
  /*!
    Copy \e other to \e this.

    \param other
      matrix to copy
   */
  inline void Update(const FixedSizeSerialDenseMatrix<rows, cols>& other);

  /// Scaled copy: \e this = \e scalarOther * \e other
  /*!
    Copy \e scalarOther * \e other to \e this.

    \param scalarOther
      scaling factor for other
    \param other
      matrix to read from
   */
  inline void Update(const double scalarOther, const FixedSizeSerialDenseMatrix<rows, cols>& other);

  /// Add: \e this = \e scalarThis * \e this + \e scalarOther * \e other
  /*!
    Scale by \e scalarThis and add \e scalarOther * \e other.

    \param scalarOther
      scaling factor for other
    \param other
      matrix to add
    \param scalarThis
      scaling factor for \e this
   */
  inline void Update(const double scalarOther, const FixedSizeSerialDenseMatrix<rows, cols>& other,
                     const double scalarThis);

  /// Add: \e this = \e left + \e right
  /*!
    Store \e left + \e right in this matrix.

    \param left
      first matrix to add
    \param right
      second matrix to add
   */
  inline void Update(const FixedSizeSerialDenseMatrix<rows, cols>& left,
                     const FixedSizeSerialDenseMatrix<rows, cols>& right);

  /// Add: \e this = \e scalarLeft * \e left + \e scalarRight * \e right
  /*!
    Store \e scalarLeft * \e left + \e scalarRight * \e right in \e this.

    \param scalarLeft
      scaling factor for \e left
    \param left
      first matrix to add
    \param scalarRight
      scaling factor for \e right
    \param right
      second matrix to add
   */
  inline void Update(const double scalarLeft, const FixedSizeSerialDenseMatrix<rows, cols>& left,
                     const double scalarRight, const FixedSizeSerialDenseMatrix<rows, cols>& right);

  /// Add: \e this = \e scalarThis * \e this + \e scalarLeft * \e left + \e scalarRight * \e right
  /*!
    Scale by \e scalarThis and add \e scalarLeft * \e left + \e scalarRight * \e right.

    \param scalarLeft
      scaling factor for \e left
    \param left
      first matrix to add
    \param scalarRight
      scaling factor for \e right
    \param right
      second matrix to add
    \param scalarThis
      scaling factor for \e this
   */
  inline void Update(const double scalarLeft, const FixedSizeSerialDenseMatrix<rows, cols>& left,
                     const double scalarRight, const FixedSizeSerialDenseMatrix<rows, cols>& right,
                     const double scalarThis);

  /// Calculate 1-norm
  /*!
    \return 1-norm
   */
  inline double Norm1() const;

  /// Calculate 2-norm (Euclidean norm)
  /*!
    \return 2-norm
   */
  inline double Norm2() const;

  /// Calculate inf-norm
  /*!
    \return inf-norm
   */
  inline double NormInf() const;

  /// Calculate minimum value
  /*!
    \return minimum value
   */
  inline double MinValue() const;

  /// Calculate maximum value
  /*!
    \return maximum value
   */
  inline double MaxValue() const;

  /// Calculate mean value
  /*!
    \return mean value
   */
  inline double MeanValue() const;

  /// Multiply: \e this = \e left*right
  /*!
    This is equivalent to MultiplyNN(\e left,\e right).

    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void Multiply(const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right);

  /// Multiply: \e this = \e left*right
  /*!
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyNN(const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right);

  /// Multiply: \e this = \e left*right^T
  /*!
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyNT(const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right);

  /// Multiply: \e this = \e left^T*right
  /*!
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyTN(const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right);

  /// Multiply: \e this = \e left^T*right^T
  /*!
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyTT(const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right);


  /// Multiply: \e this = \e scalarOthers * \e left*right
  /*!
    \param scalarOthers
      scalar factor for \e left*right
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void Multiply(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right);

  /// Multiply: \e this = \e scalarOthers * \e left*right
  /*!
    This is equivalent to MultiplyNN(\e scalarOthers,\e left,\e right).

    \param scalarOthers
      scalar factor for \e left*right
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyNN(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right);

  /// Multiply: \e this = \e scalarOthers * \e left*right^T
  /*!
    \param scalarOthers
      scalar factor for \e left*right^T
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyNT(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right);

  /// Multiply: \e this = \e scalarOthers * \e left^T*right
  /*!
    \param scalarOthers
      scalar factor for \e left^T*right
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyTN(const double scalarOthers, const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right);

  /// Multiply: \e this = \e scalarOthers * \e left^T*right^T
  /*!
    \param scalarOthers
      scalar factor for \e left^T*right^T
    \param left
      first factor
    \param right
      second factor
   */
  template <unsigned int inner>
  inline void MultiplyTT(const double scalarOthers, const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right);

  /// Multiply: \e this = \e scalarThis * \e this + \e scalarOthers * \e left*right
  /*!
    This is equivalent to MultiplyNN(\e scalarOthers,\e left,\e right,\e scalarThis).

    \param scalarOthers
      scalar factor for \e left*right
    \param left
      first factor
    \param right
      second factor
    \param scalarThis
      scalar factor for \e this
   */
  template <unsigned int inner>
  inline void Multiply(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right, const double scalarThis);

  /// Multiply: \e this = \e scalarThis * \e this + \e scalarOthers * \e left*right
  /*!
    \param scalarOthers
      scalar factor for \e left*right
    \param left
      first factor
    \param right
      second factor
    \param scalarThis
      scalar factor for \e this
   */
  template <unsigned int inner>
  inline void MultiplyNN(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right, const double scalarThis);

  /// Multiply: \e this = \e scalarThis * \e this + \e scalarOthers * \e left*right^T
  /*!
    \param scalarOthers
      scalar factor for \e left*right^T
    \param left
      first factor
    \param right
      second factor
    \param scalarThis
      scalar factor for \e this
   */
  template <unsigned int inner>
  inline void MultiplyNT(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right, const double scalarThis);

  /// Multiply: \e this = \e scalarThis * \e this + \e scalarOthers * \e left^T*right
  /*!
    \param scalarOthers
      scalar factor for \e left^T*right
    \param left
      first factor
    \param right
      second factor
    \param scalarThis
      scalar factor for \e this
   */
  template <unsigned int inner>
  inline void MultiplyTN(const double scalarOthers, const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right, const double scalarThis);

  /// Multiply: \e this = \e scalarThis * \e this + \e scalarOthers * \e left^T*right^T
  /*!
    \param scalarOthers
      scalar factor for \e left^T*right^T
    \param left
      first factor
    \param right
      second factor
    \param scalarThis
      scalar factor for \e this
   */
  template <unsigned int inner>
  inline void MultiplyTT(const double scalarOthers, const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right, const double scalarThis);

  /// Write \e this to \e out
  /*!
    Write a readable representation of \e this to \e out. This function is called by
    \e out << *\e this.

    \param out
      out stream
   */
  void Print(std::ostream& out) const;

  /// = operator
  /*!
    Copy data from \e other to \e this, equivalent to Update(other).
    \param other
      matrix to get data from
   */
  inline FixedSizeSerialDenseMatrix<rows,cols> operator=(const FixedSizeSerialDenseMatrix<rows,cols>& other);

  /// = operator for double
  /*!
    Fill with double \e other, same as PutScalar(other).

    \param other
      scalar value
   */
  inline FixedSizeSerialDenseMatrix<rows,cols> operator=(const double other);

  /// == operator
  /*!
    Compare \e this with \e other.

    \param other
      scalar value
   */
  inline bool operator==(const FixedSizeSerialDenseMatrix<rows,cols>& other) const;

  /// != operator
  /*!
    Compare \e this with \e other.

    \param other
      scalar value
   */
  inline bool operator!=(const FixedSizeSerialDenseMatrix<rows,cols>& other) const;

  /// Access data
  /*!
    Return value in row \e r and column \e c.

    \param r
      row index
    \param c
      column index
   */
  inline double& operator()(unsigned int r, unsigned int c);

  /// Access data
  /*!
    Return value in row \e r and column \e c.

    \param r
      row index
    \param c
      column index
   */
  inline const double& operator()(unsigned int r, unsigned int c) const;

  /// Access data
  /*!
    Return value in row \e r. This works only for Matrices with cols==1 or rows==1 (vectors),
    otherwise a compile time error is raised.

    \param r
      index
   */
  inline double& operator()(unsigned int r);  // for vectors, with check at compile-time

  /// Access data
  /*!
    Return value in row \e r. This works only for Matrices with cols==1 or rows==1 (vectors),
    otherwise a compile time error is raised.

    \param r
      index
   */
  inline const double& operator()(unsigned int r) const;
};

template <unsigned int cols, unsigned int rows>
std::ostream& operator<<(std::ostream& out, const FixedSizeSerialDenseMatrix<rows,cols>& matrix);

// Constructors

template<unsigned int rows, unsigned int cols>
FixedSizeSerialDenseMatrix<rows, cols>::FixedSizeSerialDenseMatrix(bool setzero) : data_(0), isview_(false)
{
  DENSEERROR::Checker<(rows != 0) and (cols != 0)>::Matrix_dimensions_cannot_be_zero();
  data_ = new double[rows*cols];
  if (setzero)
    DENSEFUNCTIONS::clearMatrix<rows,cols>(data_);
}

template<unsigned int rows, unsigned int cols>
FixedSizeSerialDenseMatrix<rows, cols>::FixedSizeSerialDenseMatrix(double* d, bool view) :
    data_(0), isview_(view)
{
  DENSEERROR::Checker<(rows != 0) and cols != 0>::Matrix_dimensions_cannot_be_zero();
  if (isview_) {
    data_ = d;
  }
  else {
    data_ = new double[rows*cols];
    std::memcpy(data_, d, rows*cols*sizeof(double));
  }
}

template<unsigned int rows, unsigned int cols>
FixedSizeSerialDenseMatrix<rows, cols>::FixedSizeSerialDenseMatrix(FixedSizeSerialDenseMatrix<rows,cols>& source, bool view) :
    data_(0), isview_(view)
{
  DENSEERROR::Checker<(rows != 0) and cols != 0>::Matrix_dimensions_cannot_be_zero();
  if (isview_) {
    data_ = source.data_;
  }
  else {
    data_ = new double[rows*cols];
    std::memcpy(data_, source.data_, rows*cols*sizeof(double));
  }
}

template<unsigned int rows, unsigned int cols>
FixedSizeSerialDenseMatrix<rows, cols>::FixedSizeSerialDenseMatrix(const FixedSizeSerialDenseMatrix<rows,cols>& source) :
  data_(0), isview_(false)
{
  data_ = new double[rows*cols];
  std::memcpy(data_, source.data_, rows*cols*sizeof(double));
}

// Destructor
template<unsigned int rows, unsigned int cols>
FixedSizeSerialDenseMatrix<rows, cols>::~FixedSizeSerialDenseMatrix()
{
  if (not isview_)
    delete[] data_;
}

template<unsigned int rows, unsigned int cols>
void FixedSizeSerialDenseMatrix<rows, cols>::SetView(double* data)
{
  if (not isview_) {
    delete[] data_;
    isview_ = true;
  }
  data_ = data;
}

template<unsigned int rows, unsigned int cols>
void FixedSizeSerialDenseMatrix<rows, cols>::SetView(FixedSizeSerialDenseMatrix<rows, cols>& source)
{
  if (not isview_) {
    delete[] data_;
    isview_ = true;
  }
  data_ = source.data_;
}

template<unsigned int rows, unsigned int cols>
void FixedSizeSerialDenseMatrix<rows, cols>::SetCopy(double* data)
{
    if (isview_) {
        data_ = new double[rows*cols];
        isview_ = false;
    }
    std::memcpy(data_, data, rows*cols*sizeof(double));
}

template<unsigned int rows, unsigned int cols>
void FixedSizeSerialDenseMatrix<rows, cols>::SetCopy(FixedSizeSerialDenseMatrix<rows, cols>& source)
{
  if (isview_) {
    data_ = new double[rows*cols];
    isview_ = false;
  }
  std::memcpy(data_, source.data_, rows*cols*sizeof(double));
}

// Determinant
template<unsigned int rows, unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::Determinant() const
{
  DENSEERROR::Checker<rows == cols>::Cannot_compute_determinant_of_nonsquare_matrix();
  return DENSEFUNCTIONS::determinant<rows, cols>(data_);
}

// Invert
template<unsigned int rows, unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::Invert()
{
  DENSEERROR::Checker<rows == cols>::Cannot_compute_inverse_of_nonsquare_matrix();
  return DENSEFUNCTIONS::invert<rows, cols>(data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::Invert(const FixedSizeSerialDenseMatrix<rows,cols>& other)
{
  DENSEERROR::Checker<rows == cols>::Cannot_compute_inverse_of_nonsquare_matrix();
  return DENSEFUNCTIONS::invert<rows, cols>(data_, other.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Clear()
{
  DENSEFUNCTIONS::clearMatrix<rows,cols>(data_);
}


template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::PutScalar(const double scalar)
{
  DENSEFUNCTIONS::putScalar<rows, cols>(scalar, data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::Dot(const FixedSizeSerialDenseMatrix<rows, cols>& other)
{
    return DENSEFUNCTIONS::dot<rows, cols>(data_, other.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Abs(const FixedSizeSerialDenseMatrix<rows, cols>& other)
{
  return DENSEFUNCTIONS::abs<rows, cols>(data_, other.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Reciprocal(const FixedSizeSerialDenseMatrix<rows, cols>& other)
{
  return DENSEFUNCTIONS::reciprocal<rows, cols>(data_, other.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Scale(const double scalar)
{
    DENSEFUNCTIONS::scaleMatrix<rows, cols>(scalar, data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Update(const FixedSizeSerialDenseMatrix<rows, cols>& other)
{
  DENSEFUNCTIONS::update<rows,cols>(data_,other.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Update(const double scalarOther,
                                                          const FixedSizeSerialDenseMatrix<rows, cols>& other)
{
  DENSEFUNCTIONS::update<rows,cols>(data_,scalarOther,other.data_);
}


template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Update(const double scalarOther,
                                                            const FixedSizeSerialDenseMatrix<rows, cols>& other,
                                                            const double scalarThis)
{
    DENSEFUNCTIONS::update<rows, cols>(scalarThis, data_, scalarOther, other.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Update(const FixedSizeSerialDenseMatrix<rows, cols>& left,
                                                          const FixedSizeSerialDenseMatrix<rows, cols>& right)
{
    DENSEFUNCTIONS::update<rows, cols>(data_, left.data_, right.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Update(const double scalarLeft,
                                                          const FixedSizeSerialDenseMatrix<rows, cols>& left,
                                                          const double scalarRight,
                                                          const FixedSizeSerialDenseMatrix<rows, cols>& right)
{
    DENSEFUNCTIONS::update<rows, cols>(data_, scalarLeft, left.data_, scalarRight, right.data_);
}

template <unsigned int rows,  unsigned int cols>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Update(const double scalarLeft,
                                                            const FixedSizeSerialDenseMatrix<rows, cols>& left,
                                                            const double scalarRight,
                                                            const FixedSizeSerialDenseMatrix<rows, cols>& right,
                                                            const double scalarThis)
{
    DENSEFUNCTIONS::update<rows, cols>(scalarThis, data_, scalarLeft, left.data_, scalarRight, right.data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::Norm1() const
{
  return DENSEFUNCTIONS::norm1<rows, cols>(data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::Norm2() const
{
  return DENSEFUNCTIONS::norm2<rows, cols>(data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::NormInf() const
{
  return DENSEFUNCTIONS::normInf<rows, cols>(data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::MinValue() const
{
  return DENSEFUNCTIONS::minValue<rows, cols>(data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::MaxValue() const
{
  return DENSEFUNCTIONS::maxValue<rows, cols>(data_);
}

template <unsigned int rows,  unsigned int cols>
inline double FixedSizeSerialDenseMatrix<rows,cols>::MeanValue() const
{
  return DENSEFUNCTIONS::meanValue<rows, cols>(data_);
}


template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Multiply(const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','N'>(data_,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyNN(const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','N'>(data_,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyNT(const FixedSizeSerialDenseMatrix<rows, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','T'>(data_,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyTN(const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'T','N'>(data_,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyTT(const FixedSizeSerialDenseMatrix<inner, rows>& left,
                       const FixedSizeSerialDenseMatrix<cols, inner>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'T','T'>(data_,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
  template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Multiply(const double scalarOthers,
                                                              const FixedSizeSerialDenseMatrix<rows, inner>& left,
                                                              const FixedSizeSerialDenseMatrix<inner, cols>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','N'>(data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyNN(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<rows, inner>& left,
                                                                const FixedSizeSerialDenseMatrix<inner, cols>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','N'>(data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyNT(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<rows, inner>& left,
                                                                const FixedSizeSerialDenseMatrix<cols, inner>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','T'>(data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyTN(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<inner, rows>& left,
                                                                const FixedSizeSerialDenseMatrix<inner, cols>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'T','N'>(data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyTT(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<inner, rows>& left,
                                                                const FixedSizeSerialDenseMatrix<cols, inner>& right)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'T','T'>(data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
  template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::Multiply(const double scalarOthers,
                                                              const FixedSizeSerialDenseMatrix<rows, inner>& left,
                                                              const FixedSizeSerialDenseMatrix<inner, cols>& right,
                                                              const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','N'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyNN(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<rows, inner>& left,
                                                                const FixedSizeSerialDenseMatrix<inner, cols>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','N'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyNT(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<rows, inner>& left,
                                                                const FixedSizeSerialDenseMatrix<cols, inner>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'N','T'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyTN(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<inner, rows>& left,
                                                                const FixedSizeSerialDenseMatrix<inner, cols>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'T','N'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows, unsigned int cols>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows,cols>::MultiplyTT(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<inner, rows>& left,
                                                                const FixedSizeSerialDenseMatrix<cols, inner>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows,inner,cols,'T','T'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows,  unsigned int cols>
void FixedSizeSerialDenseMatrix<rows,cols>::Print(std::ostream& out) const
{
  out << "FixedSizeSerialDenseMatrix<" << rows << ',' << cols << '>';
  if (isview_) out << " (view to memory only)";
  if (data_ == NULL) {
    out << " data is NULL!\n";
    return;
  }
  out << "\n[";
  for (unsigned int i = 0; i < rows; ++i) {
    if (i!=0) out << ' ';
    for (unsigned int j = 0; j < cols; ++j) {
      out << data_[i+rows*j];
      if (j < cols-1) out << ", ";
    }
    if (i < rows-1) out << ",\n";
    else out << "]\n";
  }
}

/// output operator for FixedSizeSerialDenseMatrix
/*!
  Write matrix to out. This function calls matrix.Print(out).
 */
template <unsigned int cols, unsigned int rows>
std::ostream& operator<<(std::ostream& out, const FixedSizeSerialDenseMatrix<rows,cols>& matrix)
{
    matrix.Print(out);
    return out;
}


template<unsigned int rows, unsigned int cols>
inline FixedSizeSerialDenseMatrix<rows,cols>
FixedSizeSerialDenseMatrix<rows, cols>::operator=(const FixedSizeSerialDenseMatrix<rows,cols>& other)
{
  DENSEFUNCTIONS::update<rows,cols>(data_, other.data_);
  return *this;
}

template<unsigned int rows, unsigned int cols>
inline FixedSizeSerialDenseMatrix<rows,cols>
FixedSizeSerialDenseMatrix<rows, cols>::operator=(const double other)
{
  DENSEFUNCTIONS::putScalar<rows,cols>(other, data_);
}

template<unsigned int rows, unsigned int cols>
inline bool
FixedSizeSerialDenseMatrix<rows, cols>::operator==(const FixedSizeSerialDenseMatrix<rows,cols>& other) const
{
  if (data_ == other.data_) return true;
  // unfortunately memcmp does not work, because +0 and -0 are
  // different in the memory...
  for (unsigned int c = 0; c < rows*cols; ++c)
    if (data_[c] != other.data_[c]) return false;
  return true;
}

template<unsigned int rows, unsigned int cols>
inline bool
FixedSizeSerialDenseMatrix<rows, cols>::operator!=(const FixedSizeSerialDenseMatrix<rows,cols>& other) const
{
  return not (*this == other);
}

// Access operator

template<unsigned int rows, unsigned int cols>
inline double& FixedSizeSerialDenseMatrix<rows, cols>::operator()(unsigned int r, unsigned int c)
{
#ifdef DEBUG
    if (r < 0 or r >= rows or c < 0 or c>= cols)
      dserror("Indices %i,%i out of range in FixedSizeSerialDenseMatrix<%i,%i>.",r,c,rows,cols);
#endif
  return data_[r+c*rows];
}

template<unsigned int rows, unsigned int cols>
inline const double& FixedSizeSerialDenseMatrix<rows, cols>::operator()(unsigned int r, unsigned int c) const
{
#ifdef DEBUG
    if (r < 0 or r >= rows or c < 0 or c>= cols)
      dserror("Indices %i,%i out of range in FixedSizeSerialDenseMatrix<%i,%i>.",r,c,rows,cols);
#endif
  return data_[r+c*rows];
}

template<unsigned int rows, unsigned int cols>
inline double& FixedSizeSerialDenseMatrix<rows, cols>::operator()(unsigned int r)
{
  DENSEERROR::Checker<(cols==1) or (rows==1)>::Cannot_call_1D_access_function_on_2D_matrix();
#ifdef DEBUG
    if (r < 0 or r >= (cols == 1 ? rows : cols))
      dserror("Index %i out of range in FixedSizeSerialDenseMatrix<%i,%i>.",r,rows,cols);
#endif
  return data_[r];
}

template<unsigned int rows, unsigned int cols>
inline const double& FixedSizeSerialDenseMatrix<rows, cols>::operator()(unsigned int r) const
{
    DENSEERROR::Checker<(cols==1) or (rows==1)>::Cannot_call_1D_access_function_on_2D_matrix();
#ifdef DEBUG
    if (r < 0 or r >= (cols == 1 ? rows : cols)-1)
      dserror("Index %i out of range in FixedSizeSerialDenseMatrix<%i,%i>.",r,rows,cols);
#endif
    return data_[r];
}

/// Solver for fixed size serial dense matrices
/*!
  This solver is intended to provide the funcionality of
  Epetra_SerialDenseSolver for fixed size matrices. So far only a
  subset (only the equilibration and transpose flags are available) is
  implemented for it is all that was needed. All the code of this
  solver is directly based on the Epetra solver, but with the attempt
  to simplify it and to avoid invalid states. This simplification
  might make it necessary to rewrite the class once more functionality
  is needed.

  The first two template argument specify the size of the matrix,
  although it is expected to be square. The third argument is the
  number of columns of the 'vectors'.

  \author Martin Kuettler
  \date 09/08
 */
template <unsigned int rows, unsigned int cols, unsigned int dim_rhs=1>
class FixedSizeSerialDenseSolver
{
private:

  // we do not need these functions
  FixedSizeSerialDenseSolver(const FixedSizeSerialDenseSolver<rows, cols,dim_rhs>&);
  FixedSizeSerialDenseSolver& operator=(const FixedSizeSerialDenseSolver<rows, cols,dim_rhs>&);

  /// wrapper for the LAPACK functions
  static Epetra_LAPACK lapack_;
  /// wrapper for the BLAS functions
  static Epetra_BLAS blas_;

  /// the matrix we got
  FixedSizeSerialDenseMatrix<rows,cols>* matrix_;
  /// the vector of unknowns
  FixedSizeSerialDenseMatrix<cols,dim_rhs>* vec_X_;
  /// the right hand side vector
  FixedSizeSerialDenseMatrix<rows,dim_rhs>* vec_B_;

  /// some storage for LAPACK
  std::vector<int> pivot_vec_;
  /// vector used for equilibration
  std::vector<double> r_;
  /// vector used for equilibration
  std::vector<double> c_;

  /// do we want to equilibrate?
  bool equilibrate_;
  /// should the matrix be used transposed?
  bool transpose_;
  /// is the matrix factored?
  bool factored_;
  /// is the matrix inverted?
  bool inverted_;
  /// is the system solved?
  bool solved_;


  /// Compute equilibrate scaling
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is often, but not always, a LAPACK error code.
   */
  int ComputeEquilibrateScaling();

  /// Equilibrate matrix
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is often, but not always, a LAPACK error code.
   */
  int EquilibrateMatrix();

  /// Equilibrate right hand side vector
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is often, but not always, a LAPACK error code.
   */
  int EquilibrateRHS();

  /// Unequilibrate vector of unknowns
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is often, but not always, a LAPACK error code.
   */
  int UnequilibrateLHS();
public:
  /// Constructor
  FixedSizeSerialDenseSolver();
  /// Destructor
  ~FixedSizeSerialDenseSolver();

  /// Is matrix factored?
  /*!
    \return true if matrix is factored, false otherwise
   */
  bool IsFactored() { return factored_; }

  /// Is matrix inverted?
  /*!
    \return true if matrix is inverted, false otherwise
   */
  bool IsInverted() { return inverted_; }

  /// Is system solved?
  /*!
    \return true if system is solved, false otherwise
   */
  bool IsSolved() { return solved_; }

  /// Set the matrix
  /*!
    Set the matrix to mat.

    \param mat
      new matrix
   */
  void SetMatrix(FixedSizeSerialDenseMatrix<rows,cols>& mat);

  /// Set the vectors
  /*!
    Set the vectors, the new equation is matrix*X=B.

    \param X
      vector of unknowns
    \param B
      right hand side vector
   */
  void SetVectors(FixedSizeSerialDenseMatrix<cols,dim_rhs>& X, FixedSizeSerialDenseMatrix<rows,dim_rhs>& B);

  /// Set the equilibration
  /*!
    Set whether equilibration should be used.

    \param b
      new value for equilibrate_
   */
  void FactorWithEquilibration(bool b);

  /// Set transpose
  /*!
    Set whether the matrix should be used tranposed.

    \param b
      new value for transpose_
   */
  void SolveWithTranspose(bool b) { transpose_ = b;}

  /// Factor the matrix
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is a LAPACK error code.
   */
  int Factor();

  /// Solve the system
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is a LAPACK error code or -100, indicating that
    the two vectors are the same, but may not be (when the matrix is
    inverted before the call to Solve).
   */
  int Solve();

  /// Invert the matrix
  /*
    \return integer error code. 0 if successful, negative
    otherwise. This is a LAPACK error code.
   */
  int Invert();
};

template<unsigned int rows, unsigned int cols, unsigned int dim_rhs>
FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::FixedSizeSerialDenseSolver() :
  matrix_(NULL),
  vec_X_(NULL),
  vec_B_(NULL),
  pivot_vec_(),
  r_(),
  c_(),
  equilibrate_(false),
  transpose_(false),
  factored_(false),
  inverted_(false),
  solved_(false)
{
  DENSEERROR::Checker<rows==cols>::Matrix_size_in_solver_must_be_square();
}

template<unsigned int rows, unsigned int cols, unsigned int dim_rhs>
FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::~FixedSizeSerialDenseSolver()
{
}

template<unsigned int rows, unsigned int cols, unsigned int dim_rhs>
void FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::SetMatrix(FixedSizeSerialDenseMatrix<rows,cols>& mat)
{
  c_.clear();
  r_.clear();
  pivot_vec_.clear();
  inverted_ = factored_ = solved_ = false;
  //vec_B_ = vec_X_ = NULL;
  matrix_ = &mat;
}

template<unsigned int rows, unsigned int cols, unsigned int dim_rhs>
void FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::SetVectors(FixedSizeSerialDenseMatrix<cols,dim_rhs>& X,
                                                               FixedSizeSerialDenseMatrix<rows,dim_rhs>& B)
{
  solved_ = false;
  vec_X_ = &X;
  vec_B_ = &B;
}

template<unsigned int rows, unsigned int cols, unsigned int dim_rhs>
void FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::FactorWithEquilibration(bool b)
{
#ifdef DEBUG
  if (factored_ or inverted_)
    dserror("Cannot set equilibration after changing the matrix with Factor() or Invert().");
#endif
  equilibrate_ = b;
}

template<unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::Factor()
{
#ifdef DEBUG
  if (inverted_)
    dserror("Cannot factor the inverted matrix.");
#endif
  if (factored_) return 0;
  int errnum = 0;
  if (equilibrate_) errnum = EquilibrateMatrix();
  if (errnum != 0) return errnum;
  if (pivot_vec_.empty()) pivot_vec_.resize(rows);
  lapack_.GETRF(rows, cols, matrix_->A(), rows, &(pivot_vec_[0]), &errnum);
  if (errnum != 0) return errnum;

  factored_ = true;
  return 0;
}

template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::Solve()
{
  int errnum = 0;
  if (equilibrate_) {
    errnum = EquilibrateRHS();
  }
  if (errnum != 0) return errnum;
#ifdef DEBUG
  if (not vec_B_ or not vec_X_)
    dserror("Both vectors must be set to solve.");
#endif

  if (inverted_) {
    if (vec_B_ == vec_X_) return -100;

    blas_.GEMM(transpose_ ? 'T':'N', 'N', cols, dim_rhs, cols, 1.0, matrix_->A(), rows, vec_B_->A(), rows, 0.0,
               vec_X_->A(), cols);
    solved_ = true;
  }
  else {
    if (!factored_) {
      errnum = Factor();
      if (errnum != 0) return errnum;
    }

    if (vec_B_ != vec_X_) *vec_X_ = *vec_B_;
    lapack_.GETRS(transpose_ ? 'T':'N', cols, dim_rhs, matrix_->A(), rows, &(pivot_vec_[0]), vec_X_->A(), cols, &errnum);
    if (errnum != 0) return errnum;
    solved_ = true;
  }
  if (equilibrate_) errnum = UnequilibrateLHS();
  if (errnum != 0) return errnum;
  return 0;
}

template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::ComputeEquilibrateScaling()
{
  if (!r_.empty()) return 0; //we already did that
  int errnum;
  double rowcnd, colcnd, amax;
  r_.resize(rows);
  c_.resize(cols);
  lapack_.GEEQU(rows, cols, matrix_->A(), rows, &(r_[0]), &(c_[0]), &rowcnd, &colcnd, &amax, &errnum);
  if (errnum != 0) return errnum;

  return 0;
}

template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::EquilibrateMatrix()
{
  int errnum = 0;
  if (r_.empty()) errnum = ComputeEquilibrateScaling();
  if (errnum != 0) return errnum;
  double* ptr = matrix_->A();
  double s1;
  for (unsigned j = 0; j < cols; ++j) {
    s1 = c_[j];
    for (unsigned i = 0; i < rows; ++i) {
      *ptr *= s1*r_[i];
      ++ptr;
    }
  }
  return 0;
}

template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::EquilibrateRHS()
{
  int errnum = 0;
  if (r_.empty()) errnum = ComputeEquilibrateScaling();
  if (errnum != 0) return errnum;
  std::vector<double>& r = transpose_ ? c_ : r_;
  double* ptr = vec_B_->A();
  for (unsigned j = 0; j < dim_rhs; ++j) {
    for (unsigned i = 0; i < cols; ++i) {
      *ptr *= r[i];
      ++ptr;
    }
  }

  return 0;
}

template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::UnequilibrateLHS()
{
  std::vector<double>& c = transpose_ ? r_ : c_;
  double* ptr = vec_X_->A();
  for (unsigned j = 0; j < dim_rhs; ++j) {
    for (unsigned i = 0; i < rows; ++i) {
      *ptr *= c[i];
      ++ptr;
    }
  }

  return 0;
}

template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
int FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::Invert()
{
  int errnum;
  if (not factored_) {
    errnum  = Factor();
    if (errnum != 0) return errnum;
  }

  int lwork = 4*cols;
  std::vector<double> work(lwork);
  lapack_.GETRI(cols, matrix_->A(), rows, &(pivot_vec_[0]), &(work[0]), &lwork, &errnum);
  if (errnum != 0) return errnum;
  inverted_ = true;
  factored_ = false;

  return 0;
}

// Initialize the static objects.
template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
Epetra_LAPACK FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::lapack_;
template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
Epetra_BLAS FixedSizeSerialDenseSolver<rows,cols,dim_rhs>::blas_;


}

#endif // FIXED_SIZE_SERIAL_DENSE_MATRIX_H
