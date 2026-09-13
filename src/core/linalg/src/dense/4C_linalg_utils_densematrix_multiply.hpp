// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_MULTIPLY_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_MULTIPLY_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_exceptions.hpp"

#include <cctype>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg::Internal
{
  /*!
   \brief Utility function to get type string of a matrix/vector object
   */
  template <typename VectorOrMatrix>
  inline std::string get_matrix_or_vector_string();

  template <>
  inline std::string get_matrix_or_vector_string<Core::LinAlg::SerialDenseMatrix>()
  {
    return "Matrix";
  }

  template <>
  inline std::string get_matrix_or_vector_string<Core::LinAlg::SerialDenseVector>()
  {
    return "Vector";
  }

  /*!
   \brief Utility function to get the correct case string of a matrix/vector object
   */
  template <typename VectorOrMatrix>
  inline std::string get_matrix_or_vector_case(char ch);

  template <>
  inline std::string get_matrix_or_vector_case<Core::LinAlg::SerialDenseMatrix>(char ch)
  {
    char c = static_cast<char>(std::toupper(ch));
    std::string s;
    s = c;
    return s;
  }

  template <>
  inline std::string get_matrix_or_vector_case<Core::LinAlg::SerialDenseVector>(char ch)
  {
    char c = static_cast<char>(std::tolower(ch));
    std::string s;
    s = c;
    return s;
  }

  template <typename VectorOrMatrix>
  inline int get_num_rows(const VectorOrMatrix& obj)
  {
    return obj.numRows();
  }

  template <typename VectorOrMatrix>
  inline int get_num_cols(const VectorOrMatrix& obj)
  {
    return obj.numCols();
  }

  template <>
  inline int get_num_rows<Core::LinAlg::SerialDenseVector>(
      const Core::LinAlg::SerialDenseVector& obj)
  {
    return obj.length();
  }

  template <>
  inline int get_num_cols<Core::LinAlg::SerialDenseVector>(
      const Core::LinAlg::SerialDenseVector& obj)
  {
    static_cast<void>(obj);
    return 1;
  }

  /*!
   \brief Utility function to get type string of transposition operation
   */
  template <bool transpose>
  inline std::string get_transpose_string();

  template <>
  inline std::string get_transpose_string<false>()
  {
    return "";
  }

  template <>
  inline std::string get_transpose_string<true>()
  {
    return "^T";
  }

  /*!
   \brief Utility function to check for error code of LINALG multiplication
   */
  template <bool transpose1, bool transpose2, typename VectorOrMatrix1, typename VectorOrMatrix2,
      typename VectorOrMatrix3>
  inline void check_error_code_in_debug(
      int errorCode, const VectorOrMatrix1& a, const VectorOrMatrix2& b, const VectorOrMatrix3& c)
  {
    FOUR_C_ASSERT(errorCode == 0,
        "Error code ({}) is returned. Something went wrong with the {}-{} multiplication {} = {}{} "
        "* {}{}. "
        "Dimensions of {}, {}, and {} are ({}x{}), ({}x{}), and ({}x{}) respectively.",
        errorCode, get_matrix_or_vector_string<VectorOrMatrix1>(),
        get_matrix_or_vector_string<VectorOrMatrix2>(),
        get_matrix_or_vector_case<VectorOrMatrix3>('c'),
        get_matrix_or_vector_case<VectorOrMatrix1>('a'), get_transpose_string<transpose1>(),
        get_matrix_or_vector_case<VectorOrMatrix2>('b'), get_transpose_string<transpose2>(),
        get_matrix_or_vector_case<VectorOrMatrix1>('a'),
        get_matrix_or_vector_case<VectorOrMatrix2>('b'),
        get_matrix_or_vector_case<VectorOrMatrix3>('c'), get_num_rows(a), get_num_cols(a),
        get_num_rows(b), get_num_cols(b), get_num_rows(c), get_num_cols(c));
  }
}  // namespace Core::LinAlg::Internal

namespace Core::LinAlg
{
  /*!
   \brief Matrix-vector multiplication c = A*b

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply(SerialDenseVector& c, const SerialDenseMatrix& A, const SerialDenseVector& b)
  {
    const int err = c.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A.base(), b.base(), 0.0);
    Internal::check_error_code_in_debug<false, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication c = alpha*A*b + beta*c

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply(double beta, SerialDenseVector& c, double alpha, const SerialDenseMatrix& A,
      const SerialDenseVector& b)
  {
    const int err =
        c.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A.base(), b.base(), beta);
    Internal::check_error_code_in_debug<false, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication C = A*b with matrix-shaped output

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param C (out):       Matrix C
   */
  inline int multiply(SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseVector& b)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A.base(), b.base(), 0.0);
    Internal::check_error_code_in_debug<false, false>(err, A, b, C);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication C = alpha*A*b + beta*C with matrix-shaped output

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param C (out):       Matrix C
   */
  inline int multiply(double beta, SerialDenseMatrix& C, double alpha, const SerialDenseMatrix& A,
      const SerialDenseVector& b)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A.base(), b.base(), beta);
    Internal::check_error_code_in_debug<false, false>(err, A, b, C);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication c = A^T*b

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply_tn(
      SerialDenseVector& c, const SerialDenseMatrix& A, const SerialDenseVector& b)
  {
    const int err = c.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A.base(), b.base(), 0.0);
    Internal::check_error_code_in_debug<true, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication c = alpha*A^T*b + beta*c

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply_tn(double beta, SerialDenseVector& c, double alpha,
      const SerialDenseMatrix& A, const SerialDenseVector& b)
  {
    const int err = c.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, A.base(), b.base(), beta);
    Internal::check_error_code_in_debug<true, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication C = A^T*b with matrix-shaped output

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param C (out):       Matrix C
   */
  inline int multiply_tn(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseVector& b)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A.base(), b.base(), 0.0);
    Internal::check_error_code_in_debug<true, false>(err, A, b, C);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication C = alpha*A^T*b + beta*C with matrix-shaped output

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param C (out):       Matrix C
   */
  inline int multiply_tn(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseVector& b)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, A.base(), b.base(), beta);
    Internal::check_error_code_in_debug<true, false>(err, A, b, C);
    return err;
  }

  /*!
   \brief Mixed multiplication C = a^T*B

   \param a (in):        Vector a
   \param B (in):        Matrix B
   \param C (out):       Matrix C
   */
  inline int multiply_tn(
      SerialDenseMatrix& C, const SerialDenseVector& a, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, a.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<true, false>(err, a, B, C);
    return err;
  }

  /*!
   \brief Mixed multiplication C = alpha*a^T*B + beta*C

   \param a (in):        Vector a
   \param B (in):        Matrix B
   \param C (out):       Matrix C
   */
  inline int multiply_tn(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseVector& a, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, a.base(), B.base(), beta);
    Internal::check_error_code_in_debug<true, false>(err, a, B, C);
    return err;
  }

  /*!
   \brief Mixed multiplication C = a^T*B^T

   \param a (in):        Vector a
   \param B (in):        Matrix B
   \param C (out):       Matrix C
   */
  inline int multiply_tt(
      SerialDenseMatrix& C, const SerialDenseVector& a, const SerialDenseMatrix& B)
  {
    const int err = C.base().multiply(Teuchos::TRANS, Teuchos::TRANS, 1.0, a.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<true, true>(err, a, B, C);
    return err;
  }

  /*!
   \brief Mixed multiplication C = alpha*a^T*B^T + beta*C

   \param a (in):        Vector a
   \param B (in):        Matrix B
   \param C (out):       Matrix C
   */
  inline int multiply_tt(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseVector& a, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::TRANS, alpha, a.base(), B.base(), beta);
    Internal::check_error_code_in_debug<true, true>(err, a, B, C);
    return err;
  }

  /*!
   \brief Outer product C = a*b^T

   \param a (in):        Vector a
   \param b (in):        Vector b
   \param C (out):       Matrix C
   */
  inline int multiply_nt(
      SerialDenseMatrix& C, const SerialDenseVector& a, const SerialDenseVector& b)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, a.base(), b.base(), 0.0);
    Internal::check_error_code_in_debug<false, true>(err, a, b, C);
    return err;
  }

  /*!
   \brief Outer product C = alpha*a*b^T + beta*C

   \param a (in):        Vector a
   \param b (in):        Vector b
   \param C (out):       Matrix C
   */
  inline int multiply_nt(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseVector& a, const SerialDenseVector& b)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::TRANS, alpha, a.base(), b.base(), beta);
    Internal::check_error_code_in_debug<false, true>(err, a, b, C);
    return err;
  }

  /*!
    \brief Matrix-matrix multiplication C = A*B

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply(SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<false, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A*B + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply(double beta, SerialDenseMatrix& C, double alpha, const SerialDenseMatrix& A,
      const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<false, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A^T*B

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tn(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<true, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A^T*B + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tn(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<true, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A*B^T

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_nt(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<false, true>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A*B^T + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_nt(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<false, true>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A^T*B^T

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tt(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err = C.base().multiply(Teuchos::TRANS, Teuchos::TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<true, true>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A^T*B^T + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tt(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<true, true>(err, A, B, C);
    return err;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
