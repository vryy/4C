// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_FIXEDSIZEMATRIX_SOLVER_HPP
#define FOUR_C_LINALG_FIXEDSIZEMATRIX_SOLVER_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_BLAS.hpp>
#include <Teuchos_LAPACK.hpp>

FOUR_C_NAMESPACE_OPEN


namespace Core::LinAlg
{
  /*!
   * This solver is intended to provide the functionality of
   * Epetra_SerialDenseSolver for fixed size matrices. So far only a
   * subset (only the equilibration and transpose flags are available) is
   * implemented for it is all that was needed. All the code of this
   * solver is directly based on the Epetra solver, but with the attempt
   * to simplify it and to avoid invalid states. This simplification
   * might make it necessary to rewrite the class once more functionality
   * is needed.
   *
   * The first two template argument specify the size of the matrix,
   * although it is expected to be square. The third argument is the
   * number of columns of the 'vectors'.
   */
  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs = 1>
  class FixedSizeSerialDenseSolver
  {
   private:
    static_assert(rows > 0, "FixedSizeSerialDenseSolver needs at least one row");
    static_assert(cols > 0, "FixedSizeSerialDenseSolver needs at least one row");
    static_assert(rows == cols, "FixedSizeSerialDenseSolver only works for square matrices");

    // we do not need these functions
    FixedSizeSerialDenseSolver(const FixedSizeSerialDenseSolver<rows, cols, dim_rhs>&);
    FixedSizeSerialDenseSolver& operator=(const FixedSizeSerialDenseSolver<rows, cols, dim_rhs>&);

    /// wrapper for the LAPACK functions
    static Teuchos::LAPACK<int, double> lapack_;
    /// wrapper for the BLAS functions
    static Teuchos::BLAS<unsigned int, double> blas_;

    /// the matrix we got
    Matrix<rows, cols, double>* matrix_;
    /// the vector of unknowns
    Matrix<cols, dim_rhs, double>* vec_x_;
    /// the right hand side vector
    Matrix<rows, dim_rhs, double>* vec_b_;

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
      otherwise. This is a LAPACK error code.
     */
    int compute_equilibrate_scaling();

    /// Equilibrate matrix
    /*
      \return integer error code. 0 if successful, negative
      otherwise. This is a LAPACK error code.
     */
    int equilibrate_matrix();

    /// Equilibrate right hand side vector
    /*
      \return integer error code. 0 if successful, negative
      otherwise. This is a LAPACK error code.
     */
    int equilibrate_rhs();

    /// Unequilibrate vector of unknowns
    /*
      \return integer error code. 0 if successful, negative
      otherwise. This is a LAPACK error code.
     */
    int unequilibrate_lhs();

   public:
    /// Constructor
    FixedSizeSerialDenseSolver();

    /// Is matrix factored?
    /*!
      \return true if matrix is factored, false otherwise
     */
    bool is_factored() { return factored_; }

    /// Is matrix inverted?
    /*!
      \return true if matrix is inverted, false otherwise
     */
    bool is_inverted() { return inverted_; }

    /// Is system solved?
    /*!
      \return true if system is solved, false otherwise
     */
    bool is_solved() { return solved_; }

    /// Set the matrix
    /*!
      Set the matrix to mat.

      \param mat
        new matrix
     */
    void set_matrix(Matrix<rows, cols, double>& mat);

    /// Set the vectors
    /*!
      Set the vectors, the new equation is matrix*X=B.

      \param X
        vector of unknowns
      \param B
        right hand side vector
     */
    void set_vectors(Matrix<cols, dim_rhs, double>& X, Matrix<rows, dim_rhs, double>& B);

    /// Set the equilibration
    /*!
      Set whether equilibration should be used.

      \param b
        new value for equilibrate_
     */
    void factor_with_equilibration(bool b);

    /// Set transpose
    /*!
      Set whether the matrix should be used tranposed.

      \param b
        new value for transpose_
     */
    void solve_with_transpose(bool b) { transpose_ = b; }

    /// Factor the matrix
    /*
      \return integer error code. 0 if successful, negative
      otherwise. This is a LAPACK error code.
     */
    int factor();

    /// Solve the system
    /*
      \return integer error code. 0 if successful, negative
      otherwise. This is a LAPACK error code or -100, indicating that
      the two vectors are the same, but may not be (when the matrix is
      inverted before the call to Solve).
     */
    int solve();

    /// invert the matrix
    /*
      \return integer error code. 0 if successful, negative
      otherwise. This is a LAPACK error code.
     */
    int invert();
  };

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::FixedSizeSerialDenseSolver()
      : matrix_(nullptr),
        vec_x_(nullptr),
        vec_b_(nullptr),
        pivot_vec_(),
        r_(),
        c_(),
        equilibrate_(false),
        transpose_(false),
        factored_(false),
        inverted_(false),
        solved_(false)
  {
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  void FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::set_matrix(Matrix<rows, cols, double>& mat)
  {
    c_.clear();
    r_.clear();
    pivot_vec_.clear();
    inverted_ = factored_ = solved_ = false;
    // vec_B_ = vec_X_ = nullptr;
    matrix_ = &mat;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  void FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::set_vectors(
      Matrix<cols, dim_rhs, double>& X, Matrix<rows, dim_rhs, double>& B)
  {
    solved_ = false;
    vec_x_ = &X;
    vec_b_ = &B;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  void FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::factor_with_equilibration(bool b)
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    FOUR_C_ASSERT(!factored_ && !inverted_,
        "Cannot set equilibration after changing the matrix with Factor() or invert().");
#endif
    equilibrate_ = b;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::factor()
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    FOUR_C_ASSERT(!inverted_, "Cannot factor the inverted matrix.");
#endif
    if (factored_) return 0;
    int errnum = 0;
    if (equilibrate_) errnum = equilibrate_matrix();
    if (errnum != 0) return errnum;
    if (pivot_vec_.empty()) pivot_vec_.resize(rows);
    lapack_.GETRF(rows, cols, matrix_->data(), rows, pivot_vec_.data(), &errnum);
    if (errnum != 0) return errnum;

    factored_ = true;
    return 0;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::solve()
  {
    int errnum = 0;
    if (equilibrate_)
    {
      errnum = equilibrate_rhs();
    }
    if (errnum != 0) return errnum;
#ifdef FOUR_C_ENABLE_ASSERTIONS
    FOUR_C_ASSERT(vec_b_ && vec_x_, "Both vectors must be set to solve.");
#endif

    if (inverted_)
    {
      if (vec_b_ == vec_x_) return -100;

      blas_.GEMM(transpose_ ? Teuchos::TRANS : Teuchos::NO_TRANS, Teuchos::NO_TRANS, cols, dim_rhs,
          cols, 1.0, matrix_->data(), rows, vec_b_->data(), rows, 0.0, vec_x_->data(), cols);
      solved_ = true;
    }
    else
    {
      if (!factored_)
      {
        errnum = factor();
        if (errnum != 0) return errnum;
      }

      if (vec_b_ != vec_x_) *vec_x_ = *vec_b_;
      lapack_.GETRS(transpose_ ? 'T' : 'N', cols, dim_rhs, matrix_->data(), rows, pivot_vec_.data(),
          vec_x_->data(), cols, &errnum);
      if (errnum != 0) return errnum;
      solved_ = true;
    }
    if (equilibrate_) errnum = unequilibrate_lhs();
    if (errnum != 0) return errnum;
    return 0;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::compute_equilibrate_scaling()
  {
    if (!r_.empty()) return 0;  // we already did that
    int errnum;
    double rowcnd, colcnd, amax;
    r_.resize(rows);
    c_.resize(cols);
    lapack_.GEEQU(
        rows, cols, matrix_->data(), rows, r_.data(), c_.data(), &rowcnd, &colcnd, &amax, &errnum);
    if (errnum != 0) return errnum;

    return 0;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::equilibrate_matrix()
  {
    int errnum = 0;
    if (r_.empty()) errnum = compute_equilibrate_scaling();
    if (errnum != 0) return errnum;
    double* ptr = matrix_->data();
    double s1;
    for (unsigned j = 0; j < cols; ++j)
    {
      s1 = c_[j];
      for (unsigned i = 0; i < rows; ++i)
      {
        *ptr *= s1 * r_[i];
        ++ptr;
      }
    }
    return 0;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::equilibrate_rhs()
  {
    int errnum = 0;
    if (r_.empty()) errnum = compute_equilibrate_scaling();
    if (errnum != 0) return errnum;
    std::vector<double>& r = transpose_ ? c_ : r_;
    double* ptr = vec_b_->data();
    for (unsigned j = 0; j < dim_rhs; ++j)
    {
      for (unsigned i = 0; i < cols; ++i)
      {
        *ptr *= r[i];
        ++ptr;
      }
    }

    return 0;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::unequilibrate_lhs()
  {
    std::vector<double>& c = transpose_ ? r_ : c_;
    double* ptr = vec_x_->data();
    for (unsigned j = 0; j < dim_rhs; ++j)
    {
      for (unsigned i = 0; i < rows; ++i)
      {
        *ptr *= c[i];
        ++ptr;
      }
    }

    return 0;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::invert()
  {
    int errnum;
    if (not factored_)
    {
      errnum = factor();
      if (errnum != 0) return errnum;
    }

    int lwork = 4 * cols;
    std::vector<double> work(lwork);
    lapack_.GETRI(cols, matrix_->data(), rows, pivot_vec_.data(), work.data(), lwork, &errnum);
    if (errnum != 0) return errnum;
    inverted_ = true;
    factored_ = false;

    return 0;
  }

  // Initialize the static objects.
  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  Teuchos::LAPACK<int, double> FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::lapack_;
  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  Teuchos::BLAS<unsigned int, double> FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::blas_;

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
