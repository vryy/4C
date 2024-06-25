/*----------------------------------------------------------------------*/
/*! \file

\brief a templated fixed size dense matrix

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_FIXEDSIZEMATRIX_HPP
#define FOUR_C_LINALG_FIXEDSIZEMATRIX_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_mathoperations.hpp"

#include <Teuchos_BLAS.hpp>
#include <Teuchos_LAPACK.hpp>

#include <cmath>
#include <cstring>
#include <iostream>
#include <ostream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  namespace DenseFunctions
  {
    /*
     * Declaration of the functions taking value_type*
     *
     */

    /// Multiplication: \e out = \e left*\e right
    /*!
      Multiply \e left and \e right and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c k denoting
      the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply(value_type_out* out, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e left*\e right
    /*!
      Multiply \e left and \e right and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c k denoting
      the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_nn(value_type_out* out, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e left*\e right^T
    /*!
      Multiply \e left and \e right^T and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c k denoting
      the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c k)x(\c j) so that \e
        right^T has size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_nt(value_type_out* out, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e left^T*\e right
    /*!
      Multiply \e left^T and \e right and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c k denoting
      the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param left
        pointer to the first factor, size (\c j)x(\c i) so that \e
        left^T has size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_tn(value_type_out* out, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e left^T*\e right^T
    /*!
      Multiply \e left^T and \e right^T and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c k denoting
      the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param left
        pointer to the first factor, size (\c j)x(\c i) so that \e
        left^T has size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c k)x(\c j) so that \e
        right^T has size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_tt(value_type_out* out, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e infac * \e left*\e right
    /*!
      Multiply \e left and \e right, scale the result by \e infac and store
      it in \e out. This function takes three template parameters \c
      i, \c j and \c k denoting the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right);

    /// Multiplication: \e out = \e infac * \e left*\e right
    /*!
      Multiply \e left and \e right, scale the result by \e infac and store
      it in \e out. This function takes three template parameters \c
      i, \c j and \c k denoting the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_nn(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right);

    /// Multiplication: \e out = \e infac * \e left*\e right^T
    /*!
      Multiply \e left and \e right^T, scale the result by \e infac and store
      it in \e out. This function takes three template parameters \c
      i, \c j and \c k denoting the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right^T
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c k)x(\c j) so that \e
        right^T has size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_nt(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right);

    /// Multiplication: \e out = \e infac * \e left^T*\e right
    /*!
      Multiply \e left^T and \e right, scale the result by \e infac and store
      it in \e out. This function takes three template parameters \c
      i, \c j and \c k denoting the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left^T*right
      \param left
        pointer to the first factor, size (\c j)x(\c i) so that \e
        left^T has size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_tn(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right);

    /// Multiplication: \e out = \e infac * \e left^T*\e right^T
    /*!
      Multiply \e left^T and \e right^T, scale the result by \e infac and store
      it in \e out. This function takes three template parameters \c
      i, \c j and \c k denoting the sizes of the matrices.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left^T*right^T
      \param left
        pointer to the first factor, size (\c j)x(\c i) so that \e
        left^T has size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c k)x(\c j) so that \e
        right^T has size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_tt(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right
    /*!
      Scale \e out by \e outfac and add \e left*\e right scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

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
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right
    /*!
      Scale \e out by \e outfac and add \e left*\e right scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

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
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_nn(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right^T
    /*!
      Scale \e out by \e outfac and add \e left*\e right^T scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right^T
      \param left
        pointer to the first factor, size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c k)x(\c j) so that \e
        right^T has size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_nt(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left^T*\e right
    /*!
      Scale \e out by \e outfac and add \e left^T*\e right scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left^T*right
      \param left
        pointer to the first factor, size (\c j)x(\c i) so that \e
        left^T has size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_tn(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left^T*\e right^T
    /*!
      Scale \e out by \e outfac and add \e left^T*\e right^T scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left^T*right^T
      \param left
        pointer to the first factor, size (\c j)x(\c i) so that \e
        left^T has size (\c i)x(\c j)
      \param right
        pointer to the second factor, size (\c k)x(\c j) so that \e
        right^T has size (\c j)x(\c k)
     */
    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_tt(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right);

    /// invert matrix: \e out = inv(\e in)
    /*!
      invert the matrix \e in and store the result in \e out. To keep a
      common interface there are two template parameters \c i and \c j, but
      they must be the same number. The sizes of \e in and \e out are
      expected to be (\c i)x(\c j), and they must be square.

      \param out
        pointer to the memory the result should be stored in, size (\c i)x(\c j)
      \param in
        pointer to the matrix to be inverted, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(value_type* out, const value_type* in);

    /// invert matrix: \e mat = inv(\e mat)
    /*!
      invert the matrix \e mat in place. To keep a common interface there
      are two template parameters \c i and \c j, but they must be the same
      number. The size of \e mat is expected to be (\c i)x(\c j), and it must
      be square.

      \param mat
        pointer to the matrix to be inverted in place, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(value_type* mat);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type determinant(const value_type* mat);

    /// Copy: \e out = \e in
    /*!
      Copy \e in to \e out. This function takes two template parameters \c i and
      \c j denoting the sizes of the matrices.

      \param out
        pointer to the result matrix, size (\c i)x(\c j)
      \param in
        pointer to the matrix to be copied, size (\c i)x(\c j)
     */
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_in>
    inline void update(value_type_out* out, const value_type_in* in);

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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_infac,
        class value_type_in>
    inline void update(value_type_out* out, const value_type_infac infac, const value_type_in* in);

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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_outfac,
        class value_type_infac, class value_type_in>
    inline void update(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_in* in);

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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_left,
        class value_type_right>
    inline void update(
        value_type_out* out, const value_type_left* left, const value_type_right* right);

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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_leftfac,
        class value_type_left, class value_type_rightfac, class value_type_right>
    inline void update(value_type_out* out, const value_type_leftfac leftfac,
        const value_type_left* left, const value_type_rightfac rightfac,
        const value_type_right* right);

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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_outfac,
        class value_type_leftfac, class value_type_left, class value_type_rightfac,
        class value_type_right>
    inline void update(const value_type_outfac outfac, value_type_out* out,
        const value_type_leftfac leftfac, const value_type_left* left,
        const value_type_rightfac rightfac, const value_type_right* right);

    /// Transposed copy: \e out = \e in^T
    /*!
      Copy transposed \e in to \e out. This function takes two template parameters \c i and
      \c j denoting the sizes of the matrices.

      \param out
        pointer to the result matrix, size (\c i)x(\c j)
      \param in
        pointer to the matrix to be copied, size (\c j)x(\c i)
     */
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_in>
    inline void update_t(value_type_out* out, const value_type_in* in);

    /// Scaled transposed copy: \e out = \e infac * \e in^T
    /*!
      Scale \e in by \e infac and store the transposed result in \e out. This function takes two
      template parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        pointer to the result matrix, size (\c i)x(\c j)
      \param infac
        scalar to multiply with \e in
      \param in
        pointer to the matrix to read from, size (\c j)x(\c i)
     */
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_infac,
        class value_type_in>
    inline void update_t(
        value_type_out* out, const value_type_infac infac, const value_type_in* in);

    /// Transposed addition: \e out = \e outfac * \e out + \e infac * \e in^T
    /*!
      Scale \e out by \e outfac and add \e infac * \e in^T to it. This function
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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_outfac,
        class value_type_infac, class value_type_in>
    inline void update_t(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_in* in);

    /// Multiply element-wise, \e out(m,n) = \e out(m,n)*\e in(m,n)
    /*!
      Multiply \e out and \e in, storing the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        pointer to first factor and result, size (\c i)x(\c j)
      \param in
        pointer to second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(value_type* out, const value_type* in);

    /// Multiply element-wise, \e out(m,n) = \e fac*\e out(m,n)*\e in(m,n)
    /*!
      Multiply \e out and \e in, scale by \e fac and store the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param fac
        scaling factor for the product
      \param out
        pointer to first factor and result, size (\c i)x(\c j)
      \param in
        pointer to second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type fac, value_type* out, const value_type* in);

    /// Multiply element-wise, \e out(m,n) = \e left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right and store the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        pointer to result, size (\c i)x(\c j)
      \param left
        pointer to first factor, size (\c i)x(\c j)
      \param right
        pointer to second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(
        value_type* out, const value_type* left, const value_type* right);

    /// Multiply element-wise, \e out(m,n) = \e infac*\e left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right, scale by \e infac and store the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        pointer to result, size (\c i)x(\c j)
      \param infac
        scaling factor
      \param left
        pointer to first factor, size (\c i)x(\c j)
      \param right
         pointer to second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(
        value_type* out, const value_type infac, const value_type* left, const value_type* right);

    /// Multiply element-wise, \e out(m,n) = \e outfac*\e out(m,n) + \e infac*\e left(m,n)*\e
    /// right(m,n)
    /*!
      Multiply \e left and \e right, scale by \e infac and add the result to \e out, scaled by \e
      outfac. This function takes two template parameters, unsigned ints \c i and \c j denoting the
      size of the matrices. \param outfac scaling factor for \e out \param out pointer to result,
      size (\c i)x(\c j) \param infac scaling factor the product \param left pointer to first
      factor, size (\c i)x(\c j) \param right pointer to second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type outfac, value_type* out,
        const value_type infac, const value_type* left, const value_type* right);

    /// Divide element-wise, \e out(m,n) = \e out(m,n)/\e in(m,n)
    /*!
      Devide \e out by \e in, storing the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        pointer to dividend and result, size (\c i)x(\c j)
      \param in
        pointer to divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(value_type* out, const value_type* in);

    /// Divide element-wise, \e out(m,n) = \e fac*\e out(m,n)/\e in(m,n)
    /*!
      Divide \e out by \e in, scale by \e fac and store the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param fac
        scaling factor for the product
      \param out
        pointer to dividend and result, size (\c i)x(\c j)
      \param in
        pointer to divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type fac, value_type* out, const value_type* in);

    /// Divide element-wise, \e out(m,n) = \e left(m,n)/\e right(m,n)
    /*!
      Divide \e left by \e right and store the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        pointer to result, size (\c i)x(\c j)
      \param left
        pointer to dividend, size (\c i)x(\c j)
      \param right
        pointer to divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(
        value_type* out, const value_type* left, const value_type* right);

    /// Divide element-wise, \e out(m,n) = \e infac*\e left(m,n)/\e right(m,n)
    /*!
      Divide \e left by \e right, scale by \e infac and store the result in \e out.
      This function takes two template parameters, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        pointer to result, size (\c i)x(\c j)
      \param infac
        scaling factor
      \param left
        pointer to dividend, size (\c i)x(\c j)
      \param right
         pointer to divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(
        value_type* out, const value_type infac, const value_type* left, const value_type* right);

    /// Divide element-wise, \e out(m,n) = \e outfac*\e out(m,n) + \e infac*\e left(m,n)/\e
    /// right(m,n)
    /*!
      Divide \e left by \e right, scale by \e infac and add the result to \e out, scaled by \e
      outfac. This function takes two template parameters, unsigned ints \c i and \c j denoting the
      size of the matrices. \param outfac scaling factor for \e out \param out pointer to result,
      size (\c i)x(\c j) \param infac scaling factor the product \param left pointer to dividend,
      size (\c i)x(\c j) \param right pointer to divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type outfac, value_type* out, const value_type infac,
        const value_type* left, const value_type* right);

    /// Scale matrix
    /*!
      Scale \e mat by \e fac. This function takes
      two template parameters \c i and \c j denoting the size of \e mat.

      \param fac
        scalar to multiply with \e mat
      \param mat
        pointer to the matrix, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void scale_matrix(const value_type factor, value_type* mat);

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
    template <class value_type_out, unsigned int i, unsigned int j, class value_type_left,
        class value_type_right>
    inline value_type_out dot(const value_type_left* left, const value_type_right* right);

    /// Set matrix to zero
    /*!
      Set matrix \e mat to zero. This function takes two template
      parameters i and j denoting the size of the matrix.

      This is the same as \e put_scalar<\c i, \c j>(0.0, \e mat), but it should be faster.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void clear_matrix(value_type* mat);

    /// Fill matrix with scalar value
    /*!
      Set every number in \e mat to \e scalar. This function takes two template
      parameters \c i and \c j denoting the size of the matrix.

      \param scalar
        scalar value to be set
      \param mat
        pointer to the matrix, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void put_scalar(const value_type scalar, value_type* mat);

    /// Calculate absolut values of a matrix
    /*!
      Fill \e out with the absolute values from \e in. This function takes two
      template parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        pointer to the matrix to be set, size (\c i)x(\c j)
      \param in
        pointer to the matrix the values are read from, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void abs(value_type* out, const value_type* in);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void reciprocal(value_type* out, const value_type* in);

    /// 1-norm
    /*!
      This function computes the norm of the whole matrix. It returns
      a different result than Core::LinAlg::SerialDenseMatrix::Base::OneNorm(),
      which returns the maximum of the norms of the columns.
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return 1-norm of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm1(const value_type* mat);

    /// 2-norm (Euclidean norm)
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return 2-norm of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm2(const value_type* mat);

    /// Inf-norm
    /*!
      This function does not do the same as Core::LinAlg::SerialDenseMatrix::Base::InfNorm().
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return inf-norm of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm_inf(const value_type* mat);

    /// Minimum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return minimum value of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type min_value(const value_type* mat);

    /// Maximum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return maximum value of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type max_value(const value_type* mat);

    /// Mean value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        pointer to the matrix, size (\c i)x(\c j)
      \return mean value of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type mean_value(const value_type* mat);

    /*
     * Declaration of the functions taking Core::LinAlg::SerialDenseMatrix::Base
     *
     */


    /// Multiplication: \e out = \e left*\e right
    /*!
      Multiply \e left and \e right and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c
      k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e left*\e right
    /*!
      Multiply \e left and \e right and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c
      k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nn(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e left*\e right^T
    /*!
      Multiply \e left and \e right^T and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c
      k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c k)x(\c j) so that \e right^T has size(\c
        j)x(\e k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nt(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e left^T*\e right
    /*!
      Multiply \e left^T and \e right and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c
      k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param left
        first factor, size (\c j)x(\c i) so that \e left^T has size(\c
        i)x(\e j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tn(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e left^T*\e right^T
    /*!
      Multiply \e left^T and \e right^T and store the result in \e out. This
      function takes three template parameters \c i, \c j and \c
      k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param left
        first factor, size (\c j)x(\c i) so that \e left^T has size(\c
        i)x(\e j)
      \param right
        second factor, size (\c k)x(\c j) so that \e right^T has size(\c
        j)x(\e k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tt(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e infac * \e left*\e right
    /*!
      Multiply \e left and \e right, scale the result by \e infac and store
      it in \e out. This function takes three template
      parameters \c i, \c j and \c k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e infac * \e left*\e right
    /*!
      Multiply \e left and \e right, scale the result by \e infac and store
      it in \e out. This function takes three template
      parameters \c i, \c j and \c k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nn(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e infac * \e left*\e right^T
    /*!
      Multiply \e left and \e right^T, scale the result by \e infac and store
      it in \e out. This function takes three template
      parameters \c i, \c j and \c k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c k)x(\c j) so that \e right^T has size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nt(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e infac * \e left^T*\e right
    /*!
      Multiply \e left^T and \e right, scale the result by \e infac and store
      it in \e out. This function takes three template
      parameters \c i, \c j and \c k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c j)x(\c i) so that \e left^T has size(\c
        i)x(\e j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tn(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e infac * \e left^T*\e right^T
    /*!
      Multiply \e left^T and \e right^T, scale the result by \e infac and store
      it in \e out. This function takes three template
      parameters \c i, \c j and \c k denoting the sizes of the matrices.

      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c j)x(\c i) so that \e left^T has size(\c
        i)x(\e j)
      \param right
        second factor, size (\c k)x(\c j) so that \e right^T has size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tt(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right
    /*!
      Scale \e out by \e outfac and add \e left*\e right scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

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
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right
    /*!
      Scale \e out by \e outfac and add \e left*\e right scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

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
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nn(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left*\e right^T
    /*!
      Scale \e out by \e outfac and add \e left*\e right^T scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c k)x(\c j) so that \e right^T has size
        (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nt(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left^T*\e right
    /*!
      Scale \e out by \e outfac and add \e left^T*\e right scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c j)x(\c i) so that \e left^T has size
        (\c i)x(\c j)
      \param right
        second factor, size (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tn(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiplication: \e out = \e outfac * \e out + \e infac * \e left^T*\e right^T
    /*!
      Scale \e out by \e outfac and add \e left^T*\e right^T scaled by \e
      infac. This function takes three template parameters \c i, \c j
      and \c k denoting the sizes of the matrices.

      \param outfac
        scalar to multiply with \e out
      \param out
        matrix the result should be stored in, size (\c i)x(\c k)
      \param infac
        scalar to muliply with \e left*right
      \param left
        first factor, size (\c j)x(\c i) so that \e left^T has size
        (\c i)x(\c j)
      \param right
        second factor, size (\c k)x(\c j) so that \e right^T has size
        (\c j)x(\c k)
     */
    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tt(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// invert matrix: \e out = inv(\e in)
    /*!
      invert the matrix \e in and store the result in \e out. To keep a
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
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& in);

    /// invert matrix: \e mat = inv(\e mat)
    /*!
      invert the matrix \e mat in place. To keep a common interface there
      are two template parameters \c i and \c j, but they must be the same
      number. The size of \e mat is expected to be (\c i)x(\c j), and it must
      be square.

      \note This function only works for matrices with sizes up to
      3x3. For larger matrices use the FixedSizeSerialDenseSolver.

      \param mat
        matrix to be inverted in place, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(Core::LinAlg::SerialDenseMatrix::Base& mat);

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

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type determinant(const Core::LinAlg::SerialDenseMatrix::Base& mat);

    /// Copy: \e out = \e in
    /*!
      Copy \e in to \e out. This function takes two template parameters \c i and
      \c j denoting the sizes of the matrices.

      \param out
        result matrix, size (\c i)x(\c j)
      \param in
        matrix to be copied, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& in);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& in);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void update(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& in);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type leftfac,
        const Core::LinAlg::SerialDenseMatrix::Base& left, const value_type rightfac,
        const Core::LinAlg::SerialDenseMatrix::Base& right);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void update(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type leftfac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const value_type rightfac, const Core::LinAlg::SerialDenseMatrix::Base& right);

    /// Multiply element-wise, \e out(m,n) = \e out(m,n)*\e in(m,n)
    /*!
      Multiply \e out and \e in, storing the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        first factor and result, size (\c i)x(\c j)
      \param in
        second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(
        Core::LinAlg::SerialDenseMatrix::Base out, const Core::LinAlg::SerialDenseMatrix::Base in);

    /// Multiply element-wise, \e out(m,n) = \e fac*\e out(m,n)*\e in(m,n)
    /*!
      Multiply \e out and \e in, scale by \e fac and store the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param fac
        scaling factor for the product
      \param out
        first factor and result, size (\c i)x(\c j)
      \param in
        second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type fac,
        Core::LinAlg::SerialDenseMatrix::Base out, const Core::LinAlg::SerialDenseMatrix::Base in);

    /// Multiply element-wise, \e out(m,n) = \e left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right and store the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        result, size (\c i)x(\c j)
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(Core::LinAlg::SerialDenseMatrix::Base out,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right);

    /// Multiply element-wise, \e out(m,n) = \e infac*\e left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right, scale by \e infac and store the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        result, size (\c i)x(\c j)
      \param infac
        scaling factor
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(Core::LinAlg::SerialDenseMatrix::Base out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right);

    /// Multiply element-wise, \e out(m,n) = \e outfac*out(m,n) + \e infac*\e left(m,n)*\e
    /// right(m,n)
    /*!
      Multiply \e left and \e right, scale by \e infac and add the result to \e out, scaled by \e
      outfac. This function takes two template argumens, unsigned ints \c i and \c j denoting the
      size of the matrices. \param outfac scaling factor for \e out \param out result, size (\c
      i)x(\c j) \param infac scaling factor the product \param left first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type outfac,
        Core::LinAlg::SerialDenseMatrix::Base out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right);

    /// Multiply element-wise, \e out(m,n) = \e out(m,n)*\e in(m,n)
    /*!
      Multiply \e out and \e in, storing the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        dividend and result, size (\c i)x(\c j)
      \param in
        divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(
        Core::LinAlg::SerialDenseMatrix::Base out, const Core::LinAlg::SerialDenseMatrix::Base in);

    /// Divide element-wise, \e out(m,n) = \e fac*\e out(m,n)*\e in(m,n)
    /*!
      Divide \e out and \e in, scale by \e fac and store the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param fac
        scaling factor for the product
      \param out
        dividend and result, size (\c i)x(\c j)
      \param in
        divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type fac, Core::LinAlg::SerialDenseMatrix::Base out,
        const Core::LinAlg::SerialDenseMatrix::Base in);

    /// Divide element-wise, \e out(m,n) = \e left(m,n)*\e right(m,n)
    /*!
      Divide \e left and \e right and store the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        result, size (\c i)x(\c j)
      \param left
        dividend, size (\c i)x(\c j)
      \param right
        divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(Core::LinAlg::SerialDenseMatrix::Base out,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right);

    /// Divide element-wise, \e out(m,n) = \e infac*\e left(m,n)*\e right(m,n)
    /*!
      Divide \e left and \e right, scale by \e infac and store the result in \e out.
      This function takes two template argumens, unsigned ints \c i and \c j
      denoting the size of the matrices.
      \param out
        result, size (\c i)x(\c j)
      \param infac
        scaling factor
      \param left
        dividend, size (\c i)x(\c j)
      \param right
        divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(Core::LinAlg::SerialDenseMatrix::Base out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right);

    /// Divide element-wise, \e out(m,n) = \e outfac*\e out(m,n) + \e infac*\e left(m,n)*\e
    /// right(m,n)
    /*!
      Divide \e left by \e right, scale by \e infac and add the result to \e out, scaled by \e
      outfac. This function takes two template argumens, unsigned ints \c i and \c j denoting the
      size of the matrices. \param outfac scaling factor for \e out \param out result, size (\c
      i)x(\c j) \param infac scaling factor the product \param left dividend, size (\c i)x(\c j)
      \param right
        divisor, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type outfac,
        Core::LinAlg::SerialDenseMatrix::Base out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right);

    /// Scale matrix
    /*!
      Scale \e mat by \e scalar. This function takes
      two template parameters \c i and \c j denoting the size of \e mat.

      \param scalar
        scalar to multiply with \e mat
      \param mat
        matrix, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void scale_matrix(const value_type scalar, Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      mat.scale(scalar);
    }

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
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type dot(const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
      return dot<value_type, i, j>(left.values(), right.values());
    }

    /// Set matrix to zero
    /*!
      Set matrix \e mat to zero. This function takes two template
      parameters i and j denoting the size of the matrix.

      This is the same as \e put_scalar<\c i,\c j>(0.0, \e mat), but it should be faster.

      \param mat
        matrix, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void clear_matrix(Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      mat.putScalar(0.0);
    }

    /// Fill matrix with scalar value
    /*!
      Set every number in \e mat to \e scalar. This function takes two template
      parameters \c i and \c j denoting the size of the matrix.

      \param scalar
        scalar value to be set
      \param mat
        matrix, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void put_scalar(const value_type scalar, Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      mat.putScalar(scalar);
    }

    /// Calculate absolut values of a matrix
    /*!
      Fill \e out with the absolute values from \e in. This function takes two
      template parameters \c i and \c j denoting the sizes of the matrices.

      \param out
        matrix to be set, size (\c i)x(\c j)
      \param in
        matrix the values are read from, size (\c i)x(\c j)
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline void abs(Core::LinAlg::SerialDenseMatrix::Base& dest,
        const Core::LinAlg::SerialDenseMatrix::Base& src);

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
    template <class value_type, unsigned int i, unsigned int j>
    inline void reciprocal(Core::LinAlg::SerialDenseMatrix::Base& dest,
        const Core::LinAlg::SerialDenseMatrix::Base& src);

    /// 1-norm
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return 1-norm of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm1(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      return norm1<value_type, i, j>(mat.values());
    }

    /// 2-norm (Euclidean norm)
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return 2-norm of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm2(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      return norm2<value_type, i, j>(mat.values());
    }

    /// Inf-norm
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return inf-norm of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm_inf(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      return norm_inf<value_type, i, j>(mat.values());
    }

    /// Minimum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return minimum value of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type min_value(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      return min_value<value_type, i, j>(mat.values());
    }

    /// Maximum value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return maximum value of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type max_value(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      return max_value<value_type, i, j>(mat.values());
    }

    /// Mean value of a matrix
    /*!
      The template arguments \c i and \c j are the size of the matrix.

      \param mat
        matrix, size (\c i)x(\c j)
      \return mean value of \e mat
     */
    template <class value_type, unsigned int i, unsigned int j>
    inline value_type mean_value(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
      return mean_value<value_type, i, j>(mat.values());
    }

    /*
     * Definitions of the functions taking value_type*
     *
     */

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply(
        value_type_out* out, const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3];
          }
          *out = tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_nn(
        value_type_out* out, const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3];
          }
          *out = tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_nt(
        value_type_out* out, const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < k; ++c1)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3 * k];
          }
          *out = tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_tn(
        value_type_out* out, const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i * j; c2 += j)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3] * right[c1 + c3];
          }
          *out = tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_left, class value_type_right>
    inline void multiply_tt(
        value_type_out* out, const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < k; ++c1)
      {
        for (unsigned int c2 = 0; c2 < i * j; c2 += j)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3] * right[c1 + c3 * k];
          }
          *out = tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3];
          }
          *out = infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_nn(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3];
          }
          *out = infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_nt(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < k; ++c1)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3 * k];
          }
          *out = infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_tn(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i * j; c2 += j)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3] * right[c1 + c3];
          }
          *out = infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_infac, class value_type_left, class value_type_right>
    inline void multiply_tt(value_type_out* out, const value_type_infac infac,
        const value_type_left* const left, const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < k; ++c1)
      {
        for (unsigned int c2 = 0; c2 < i * j; c2 += j)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3] * right[c1 + c3 * k];
          }
          *out = infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3];
          }
          *out = (*out) * outfac + infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_nn(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3];
          }
          *out = (*out) * outfac + infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_nt(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < k; ++c1)
      {
        for (unsigned int c2 = 0; c2 < i; ++c2)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3 * i] * right[c1 + c3 * k];
          }
          *out = (*out) * outfac + infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_tn(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < j * k; c1 += j)
      {
        for (unsigned int c2 = 0; c2 < i * j; c2 += j)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3] * right[c1 + c3];
          }
          *out = (*out) * outfac + infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, unsigned int k,
        class value_type_outfac, class value_type_infac, class value_type_left,
        class value_type_right>
    inline void multiply_tt(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_left* const left,
        const value_type_right* const right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      for (unsigned int c1 = 0; c1 < k; ++c1)
      {
        for (unsigned int c2 = 0; c2 < i * j; c2 += j)
        {
          value_type_out tmp = left[c2] * right[c1];
          for (unsigned int c3 = 1; c3 < j; ++c3)
          {
            tmp += left[c2 + c3] * right[c1 + c3 * k];
          }
          *out = (*out) * outfac + infac * tmp;
          ++out;
        }
      }
    }

    template <class value_type>
    inline value_type invert1x1(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      const value_type det = in[0];
      if (det == 0.0) FOUR_C_THROW("determinant of 1x1 matrix is zero");
      out[0] = 1.0 / in[0];
      return det;
    }

    template <class value_type>
    inline value_type invert2x2(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      const value_type det = in[0] * in[1 + 1 * 2] - in[1] * in[1 * 2];
      if (det == 0.0) FOUR_C_THROW("determinant of 2x2 matrix is zero");
      const value_type invdet = 1.0 / det;
      out[0] = invdet * in[1 + 1 * 2];
      out[1] = -invdet * in[1];
      out[1 * 2] = -invdet * in[1 * 2];
      out[1 + 1 * 2] = invdet * in[0];
      return det;
    }


    template <class value_type>
    inline value_type invert3x3(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      out[0] = in[1 + 1 * 3] * in[2 + 2 * 3] - in[2 + 1 * 3] * in[1 + 2 * 3];
      out[1] = in[2] * in[1 + 2 * 3] - in[1] * in[2 + 2 * 3];
      out[2] = in[1] * in[2 + 1 * 3] - in[2] * in[1 + 1 * 3];
      const value_type det = in[0] * out[0] + in[1 * 3] * out[1] + in[2 * 3] * out[2];
      // const value_type det = in[0]*in[1+3*1]*in[2+3*2] +
      //                   in[0+3*1]*in[1+3*2]*in[2+3*0] +
      //                   in[0+3*2]*in[1+3*0]*in[2+3*1] -
      //                   in[0+3*2]*in[1+3*1]*in[2+3*0] -
      //                   in[0+3*0]*in[1+3*2]*in[2+3*1] -
      //                   in[0+3*1]*in[1+3*0]*in[2+3*2];
      if (det == 0.0) FOUR_C_THROW("determinant of 3x3 matrix is zero");
      const value_type invdet = 1.0 / det;
      out[0] *= invdet;
      out[1] *= invdet;
      out[2] *= invdet;
      out[1 * 3] = invdet * (in[2 + 1 * 3] * in[2 * 3] - in[1 * 3] * in[2 + 2 * 3]);
      out[1 + 1 * 3] = invdet * (in[0] * in[2 + 2 * 3] - in[2] * in[2 * 3]);
      out[2 + 1 * 3] = invdet * (in[2] * in[1 * 3] - in[0] * in[2 + 1 * 3]);
      out[2 * 3] = invdet * (in[1 * 3] * in[1 + 2 * 3] - in[1 + 1 * 3] * in[2 * 3]);
      out[1 + 2 * 3] = invdet * (in[1] * in[2 * 3] - in[0] * in[1 + 2 * 3]);
      out[2 + 2 * 3] = invdet * (in[0] * in[1 + 1 * 3] - in[1] * in[1 * 3]);
      return det;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      static_assert(i == j, "Cannot compute inverse of non-square matrix");

      switch (i)
      {
        case 1:
          return invert1x1(out, in);
        case 2:
          return invert2x2(out, in);
        case 3:
          return invert3x3(out, in);
        default:
          static_assert(i < 4, "Cannot compute inverse of matrix bigger than 3x3");
          return 0.0;
      }
    }


    template <class value_type>
    inline value_type invert1x1(value_type* mat)
    {
      const value_type det = mat[0];
      if (det == 0.0) FOUR_C_THROW("determinant of 1x1 matrix is zero");
      mat[0] = 1.0 / mat[0];
      return det;
    }

    template <class value_type>
    inline value_type invert2x2(value_type* mat)
    {
      value_type tmp;
      const value_type det = mat[0] * mat[1 + 1 * 2] - mat[1] * mat[1 * 2];
      if (det == 0.0) FOUR_C_THROW("determinant of 2x2 matrix is zero");
      const value_type invdet = 1.0 / det;
      tmp = mat[0];
      mat[0] = invdet * mat[1 + 1 * 2];
      mat[1 + 1 * 2] = invdet * tmp;
      mat[1] *= -invdet;
      mat[1 * 2] *= -invdet;
      return det;
    }


    template <class value_type>
    inline value_type invert3x3(value_type* mat)
    {
      const value_type tmp00 = mat[1 + 1 * 3] * mat[2 + 2 * 3] - mat[2 + 1 * 3] * mat[1 + 2 * 3];
      const value_type tmp10 = mat[2] * mat[1 + 2 * 3] - mat[1] * mat[2 + 2 * 3];
      const value_type tmp20 = mat[1] * mat[2 + 1 * 3] - mat[2] * mat[1 + 1 * 3];
      const value_type det = mat[0] * tmp00 + mat[1 * 3] * tmp10 + mat[2 * 3] * tmp20;
      // const value_type det = mat[0+3*0]*mat[1+3*1]*mat[2+3*2] +
      //                    mat[0+3*1]*mat[1+3*2]*mat[2+3*0] +
      //                    mat[0+3*2]*mat[1+3*0]*mat[2+3*1] -
      //                    mat[0+3*2]*mat[1+3*1]*mat[2+3*0] -
      //                    mat[0+3*0]*mat[1+3*2]*mat[2+3*1] -
      //                    mat[0+3*1]*mat[1+3*0]*mat[2+3*2];
      if (det == 0.0) FOUR_C_THROW("determinant of 3x3 matrix is zero");
      const value_type invdet = 1.0 / det;
      const value_type tmp01 = mat[1 * 3];
      const value_type tmp11 = mat[1 + 1 * 3];
      const value_type tmp12 = mat[1 + 2 * 3];
      mat[1 * 3] = invdet * (mat[2 + 1 * 3] * mat[2 * 3] - tmp01 * mat[2 + 2 * 3]);
      mat[1 + 1 * 3] = invdet * (mat[0] * mat[2 + 2 * 3] - mat[2] * mat[2 * 3]);
      mat[1 + 2 * 3] = invdet * (mat[1] * mat[2 * 3] - mat[0] * tmp12);
      mat[2 + 1 * 3] = invdet * (mat[2] * tmp01 - mat[0] * mat[2 + 1 * 3]);
      mat[2 * 3] = invdet * (tmp01 * tmp12 - tmp11 * mat[2 * 3]);
      mat[2 + 2 * 3] = invdet * (mat[0] * tmp11 - mat[1] * tmp01);
      mat[0] = invdet * tmp00;
      mat[1] = invdet * tmp10;
      mat[2] = invdet * tmp20;
      return det;
    }


    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(value_type* mat)
    {
      static_assert(i == j, "Cannot compute inverse of non-square matrix");

      switch (i)
      {
        case 1:
          return invert1x1(mat);
        case 2:
          return invert2x2(mat);
        case 3:
          return invert3x3(mat);
        default:
          static_assert(i < 4, "Cannot compute inverse of matrix bigger than 3x3");
          return 0.0;
      }
    }

    template <class value_type>
    inline value_type determinant_large_matrix(
        unsigned int i, unsigned int j, const value_type* mat)
    {
      FOUR_C_THROW("determinant_large_matrix not implemented for this value_type!");
      return 0.0;
    }

    // specialization for double as lapack routine is used
    template <>
    inline double determinant_large_matrix<double>(
        unsigned int i, unsigned int j, const double* mat)
    {
      // taken from src/linalg/linalg_utils_densematrix_eigen.cpp: Core::LinAlg::DeterminantLU,
      // only with minor changes.
      std::vector<double> tmp(i * j);
      std::copy(mat, mat + i * j, tmp.data());
      std::vector<int> ipiv(j);
      int info;

      Teuchos::LAPACK<int, double> lapack;
      lapack.GETRF(i, j, tmp.data(), i, ipiv.data(), &info);

      if (info < 0)
        FOUR_C_THROW("Lapack's dgetrf returned %d", info);
      else if (info > 0)
        return 0.0;
      double d = tmp[0];
      for (unsigned int c = 1; c < j; ++c) d *= tmp[c + i * c];
      // swapping rows of A changes the sign of the determinant, so we have to
      // undo lapack's permutation w.r.t. the determinant
      // note the fortran indexing convention in ipiv
      for (unsigned int c = 0; c < j; ++c)
        if (static_cast<unsigned>(ipiv[c]) != c + 1) d *= -1.0;
      return d;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type determinant(const value_type* mat)
    {
      static_assert(i == j, "Matrix must be square");

      switch (i)
      {
        case 1:
          return *mat;
        case 2:
          return mat[0] * mat[1 + 1 * 2] - mat[1] * mat[1 * 2];
        case 3:
          return mat[0] * (mat[1 + 1 * 3] * mat[2 + 2 * 3] - mat[2 + 1 * 3] * mat[1 + 2 * 3]) +
                 mat[1 * 3] * (mat[2] * mat[1 + 2 * 3] - mat[1] * mat[2 + 2 * 3]) +
                 mat[2 * 3] * (mat[1] * mat[2 + 1 * 3] - mat[2] * mat[1 + 1 * 3]);
        default:
          return determinant_large_matrix<value_type>(i, j, mat);
      }
    }


    /* add matrices */

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_in>
    inline void update(value_type_out* out, const value_type_in* in)
    {
      if constexpr (std::is_same_v<value_type_out, value_type_in>)
      {
#ifdef FOUR_C_DEBUG
        if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
        // std::memcpy(out, in, i*j*sizeof(value_type));
        std::copy(in, in + i * j, out);
      }
      else
      {
        update<value_type_out, i, j>(out, 1.0, in);
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_infac,
        class value_type_in>
    inline void update(value_type_out* out, const value_type_infac infac, const value_type_in* in)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_in>)
        if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out = infac * (*in);
      for (unsigned int c = 1; c < i * j; ++c) *(++out) = infac * (*(++in));
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_outfac,
        class value_type_infac, class value_type_in>
    inline void update(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_in* in)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_in>)
        if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      if (outfac > -1e-30 and outfac < 1e-30)
      {  // cannot handle this case here, because 0*nan==nan
        update<value_type_out, i, j>(out, infac, in);
        return;
      }
      *out *= outfac;
      *out += infac * (*in);
      for (unsigned int c = 1; c < i * j; ++c)
      {
        *(++out) *= outfac;
        *out += infac * (*(++in));
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_left,
        class value_type_right>
    inline void update(
        value_type_out* out, const value_type_left* left, const value_type_right* right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      *out = *left + *right;
      for (unsigned int c = 1; c < i * j; ++c) *(++out) = *(++left) + *(++right);
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_leftfac,
        class value_type_left, class value_type_rightfac, class value_type_right>
    inline void update(value_type_out* out, const value_type_leftfac leftfac,
        const value_type_left* left, const value_type_rightfac rightfac,
        const value_type_right* right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      *out = leftfac * (*left) + rightfac * (*right);
      for (unsigned int c = 1; c < i * j; ++c)
        *(++out) = leftfac * (*(++left)) + rightfac * (*(++right));
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_outfac,
        class value_type_leftfac, class value_type_left, class value_type_rightfac,
        class value_type_right>
    inline void update(const value_type_outfac outfac, value_type_out* out,
        const value_type_leftfac leftfac, const value_type_left* left,
        const value_type_rightfac rightfac, const value_type_right* right)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_left>)
        if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if constexpr (std::is_same_v<value_type_out, value_type_right>)
        if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      if (outfac > -1e-30 and outfac < 1e-30)
      {  // cannot handle this case here, because 0*nan==nan
        update<value_type_out, i, j>(out, leftfac, left, rightfac, right);
        return;
      }
      *out *= outfac;
      *out += leftfac * (*left) + rightfac * (*right);
      for (unsigned int c = 1; c < i * j; ++c)
      {
        *(++out) *= outfac;
        *out += leftfac * (*(++left)) + rightfac * (*(++right));
      }
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_in>
    inline void update_t(value_type_out* out, const value_type_in* in)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_in>)
        if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      for (unsigned int c2 = 0; c2 < j; c2 += 1)
        for (unsigned int c1 = 0; c1 < i; c1 += 1) *(out++) = in[c2 + c1 * j];
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_infac,
        class value_type_in>
    inline void update_t(value_type_out* out, const value_type_infac infac, const value_type_in* in)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_in>)
        if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      for (unsigned int c2 = 0; c2 < j; c2 += 1)
        for (unsigned int c1 = 0; c1 < i; c1 += 1) *(out++) = infac * in[c2 + c1 * j];
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_outfac,
        class value_type_infac, class value_type_in>
    inline void update_t(const value_type_outfac outfac, value_type_out* out,
        const value_type_infac infac, const value_type_in* in)
    {
#ifdef FOUR_C_DEBUG
      if constexpr (std::is_same_v<value_type_out, value_type_in>)
        if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      if (outfac > -1e-30 and outfac < 1e-30)
      {  // cannot handle this case here, because 0*nan==nan
        update_t<value_type_out, i, j>(out, infac, in);
        return;
      }
      for (unsigned int c2 = 0; c2 < j; c2 += 1)
      {
        for (unsigned int c1 = 0; c1 < i; c1 += 1)
        {
          *(out) *= outfac;
          *(out++) += infac * in[c2 + c1 * j];
        }
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out *= *in;
      for (unsigned c = 1; c < i * j; ++c) *(++out) *= *(++in);
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type fac, value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out *= fac * (*in);
      for (unsigned c = 1; c < i * j; ++c) *(++out) *= fac * (*(++in));
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(
        value_type* out, const value_type* left, const value_type* right)
    {
#ifdef FOUR_C_DEBUG
      if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      *out = (*left) * (*right);
      for (unsigned c = 1; c < i * j; ++c) *(++out) = (*(++left)) * (*(++right));
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(
        value_type* out, const value_type infac, const value_type* left, const value_type* right)
    {
#ifdef FOUR_C_DEBUG
      if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      *out = infac * (*left) * (*right);
      for (unsigned c = 1; c < i * j; ++c) *(++out) = infac * (*(++left)) * (*(++right));
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type outfac, value_type* out,
        const value_type infac, const value_type* left, const value_type* right)
    {
#ifdef FOUR_C_DEBUG
      if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      if (outfac > -1e-30 and outfac < 1e-30)
      {
        elementwise_multiply<value_type, i, j>(out, infac, left, right);
        return;
      }
      *out = outfac * (*out) + infac * (*left) * (*right);
      for (unsigned c = 1; c < i * j; ++c)
      {
        ++out;
        *out = outfac * (*out) + infac * (*(++left)) * (*(++right));
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out /= *in;
      for (unsigned c = 1; c < i * j; ++c) *(++out) /= *(++in);
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type fac, value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out = fac * (*out) / (*in);
      for (unsigned c = 1; c < i * j; ++c)
      {
        ++out;
        ++in;
        *out = fac * (*out) / (*in);
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(value_type* out, const value_type* left, const value_type* right)
    {
#ifdef FOUR_C_DEBUG
      if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      *out = (*left) / (*right);
      for (unsigned c = 1; c < i * j; ++c) *(++out) = (*(++left)) / (*(++right));
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(
        value_type* out, const value_type infac, const value_type* left, const value_type* right)
    {
#ifdef FOUR_C_DEBUG
      if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      *out = infac * (*left) / (*right);
      for (unsigned c = 1; c < i * j; ++c) *(++out) = infac * (*(++left)) / (*(++right));
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type outfac, value_type* out, const value_type infac,
        const value_type* left, const value_type* right)
    {
#ifdef FOUR_C_DEBUG
      if (out == left) FOUR_C_THROW("'out' and 'left' point to same memory location");
      if (out == right) FOUR_C_THROW("'out' and 'right' point to same memory location");
#endif
      if (outfac > -1e-30 and outfac < 1e-30)
      {
        elementwise_divide<value_type, i, j>(out, infac, left, right);
        return;
      }
      *out = outfac * (*out) + infac * (*left) / (*right);
      for (unsigned c = 1; c < i * j; ++c)
      {
        ++out;
        *out = outfac * (*out) + infac * (*(++left)) / (*(++right));
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void scale_matrix(const value_type factor, value_type* mat)
    {
      *mat *= factor;
      for (unsigned int c = 1; c < i * j; ++c) *(++mat) *= factor;
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_left,
        class value_type_right>
    inline value_type_out dot(const value_type_left* left, const value_type_right* right)
    {
      value_type_out res = (*left) * (*right);
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++left;
        ++right;
        res += (*left) * (*right);
      }
      return res;
    }

    template <class value_type_out, unsigned int i, unsigned int j, class value_type_left,
        class value_type_right>
    inline void crossproduct(
        value_type_out* out, const value_type_left* left, const value_type_right* right)
    {
#ifdef FOUR_C_DEBUG
      if (i != 3 || j != 1) FOUR_C_THROW("cross product only for 3x1 matrices available");
#endif
      out[0] = left[1] * right[2] - left[2] * right[1];
      out[1] = left[2] * right[0] - left[0] * right[2];
      out[2] = left[0] * right[1] - left[1] * right[0];
      return;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void clear_matrix(value_type* mat)
    {
      // the memset method is needed for arbitrary precision (cln) data types instead of the fill
      // method std::memset(mat,0,i*j*sizeof(value_type));
      std::fill(mat, mat + i * j, 0.0);
    }


    template <class value_type, unsigned int i, unsigned int j>
    inline void put_scalar(const value_type scalar, value_type* mat)
    {
      *mat = scalar;
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++mat;
        *mat = scalar;
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void abs(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out = *in >= 0 ? *in : -*in;
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++out;
        ++in;
        *out = *in >= 0 ? *in : -*in;
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void reciprocal(value_type* out, const value_type* in)
    {
#ifdef FOUR_C_DEBUG
      if (out == in) FOUR_C_THROW("'out' and 'in' point to same memory location");
#endif
      *out = 1.0 / (*in);
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++out;
        ++in;
        *out = 1.0 / (*in);
      }
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm1(const value_type* mat)
    {
      value_type result = *mat >= 0 ? *mat : -(*mat);
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++mat;
        result += *mat >= 0 ? *mat : -(*mat);
      }
      return result;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm2(const value_type* mat)
    {
      value_type result = (*mat) * (*mat);
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++mat;
        result += (*mat) * (*mat);
      }
      return Core::MathOperations<value_type>::sqrt(result);
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type norm_inf(const value_type* mat)
    {
      value_type result = Core::MathOperations<value_type>::abs(*mat);
      value_type tmp;
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++mat;
        tmp = Core::MathOperations<value_type>::abs(*mat);
        result = std::max(result, tmp);
      }
      return result;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type min_value(const value_type* mat)
    {
      value_type result = *mat;
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++mat;
        if (*mat < result) result = *mat;
      }
      return result;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type max_value(const value_type* mat)
    {
      value_type result = *mat;
      for (unsigned int c = 1; c < i * j; ++c)
      {
        ++mat;
        if (*mat > result) result = *mat;
      }
      return result;
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type mean_value(const value_type* mat)
    {
      value_type result = *mat;
      for (unsigned int c = 1; c < i * j; ++c) result += *(++mat);
      return result / (i * j);
    }

    /*
     * Definitions of the functions taking Core::LinAlg::SerialDenseMatrix::Base
     *
     */

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numCols() or
          left.numCols() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply<value_type, i, j, k>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nn(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numCols() or
          left.numCols() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_nn<value_type, i, j, k>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nt(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numRows() or
          left.numCols() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_nt<value_type, i, j, k>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tn(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numCols() or out.numCols() != right.numCols() or
          left.numRows() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_tn<value_type, i, j, k>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tt(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numCols() or out.numCols() != right.numRows() or
          left.numRows() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_tt<value_type, i, j, k>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numCols() or
          left.numCols() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply<value_type, i, j, k>(out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nn(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numCols() or
          left.numCols() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_nn<value_type, i, j, k>(out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nt(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numRows() or
          left.numCols() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_nt<value_type, i, j, k>(out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tn(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numCols() or out.numCols() != right.numCols() or
          left.numRows() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_tn<value_type, i, j, k>(out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tt(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numCols() or out.numCols() != right.numRows() or
          left.numRows() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_tt<value_type, i, j, k>(out.values(), infac, left.values(), right.values());
    }


    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numCols() or
          left.numCols() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply<value_type, i, j, k>(outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nn(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numCols() or
          left.numCols() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_nn<value_type, i, j, k>(outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_nt(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or out.numCols() != right.numRows() or
          left.numCols() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_nt<value_type, i, j, k>(outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tn(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numCols() or out.numCols() != right.numCols() or
          left.numRows() != right.numRows())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_tn<value_type, i, j, k>(outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j, unsigned int k>
    inline void multiply_tt(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numCols() or out.numCols() != right.numRows() or
          left.numRows() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      multiply_tt<value_type, i, j, k>(outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type invert(
        Core::LinAlg::SerialDenseMatrix::Base& out, const Core::LinAlg::SerialDenseMatrix::Base& in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != out.numCols() or in.numRows() != in.numCols() or
          out.numRows() != in.numRows())
        FOUR_C_THROW("Invalid matrix sizes for inversion, (%i,%i) = inv( (%i,%i) )", out.numRows(),
            out.numCols(), in.numRows(), in.numCols());
#endif
      return invert<value_type, i, j>(out.values(), in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline value_type determinant(const Core::LinAlg::SerialDenseMatrix::Base& mat)
    {
#ifdef FOUR_C_DEBUG
      if (mat.numRows() != mat.numCols())
        FOUR_C_THROW(
            "Invalid matrix sizes for determinant, inv( (%i,%i) )", mat.numRows(), mat.numCols());
#endif
      return determinant<value_type, i, j>(mat.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out,
        const Core::LinAlg::SerialDenseMatrix::Base& left,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or left.numRows() != right.numRows() or
          out.numCols() != left.numCols() or left.numCols() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      update<value_type, i, j>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type leftfac,
        const Core::LinAlg::SerialDenseMatrix::Base& left, const value_type rightfac,
        const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or left.numRows() != right.numRows() or
          out.numCols() != left.numCols() or left.numCols() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      update<value_type, i, j>(out.values(), leftfac, left.values(), rightfac, right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void update(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type leftfac, const Core::LinAlg::SerialDenseMatrix::Base& left,
        const value_type rightfac, const Core::LinAlg::SerialDenseMatrix::Base& right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != left.numRows() or left.numRows() != right.numRows() or
          out.numCols() != left.numCols() or left.numCols() != right.numCols())
        FOUR_C_THROW("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      update<value_type, i, j>(
          outfac, out.values(), leftfac, left.values(), rightfac, right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void update(
        Core::LinAlg::SerialDenseMatrix::Base& out, const Core::LinAlg::SerialDenseMatrix::Base& in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != in.numRows() or out.numCols() != in.numCols())
        FOUR_C_THROW("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)", out.numRows(),
            out.numCols(), in.numRows(), in.numCols());
#endif
      update<value_type, i, j>(out.values(), in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void update(Core::LinAlg::SerialDenseMatrix::Base& out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base& in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != in.numRows() or out.numCols() != in.numCols())
        FOUR_C_THROW("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)", out.numRows(),
            out.numCols(), in.numRows(), in.numCols());
#endif
      update<value_type, i, j>(out.values(), infac, in.values());
    }


    template <class value_type, unsigned int i, unsigned int j>
    inline void update(const value_type outfac, Core::LinAlg::SerialDenseMatrix::Base& out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base& in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != in.numRows() or out.numCols() != in.numCols())
        FOUR_C_THROW("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)", out.numRows(),
            out.numCols(), in.numRows(), in.numCols());
#endif
      update<value_type, i, j>(outfac, out.values(), infac, in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(
        Core::LinAlg::SerialDenseMatrix::Base out, const Core::LinAlg::SerialDenseMatrix::Base in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or in.numRows() != i or in.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_multiply<%i,%i>, (%i,%i) *= (%i,%i)", i,
            j, out.numRows(), out.numCols(), in.numRows(), in.numCols());
#endif
      elementwise_multiply<value_type, i, j>(out.values(), in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type fac,
        Core::LinAlg::SerialDenseMatrix::Base out, const Core::LinAlg::SerialDenseMatrix::Base in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or in.numRows() != i or in.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_multiply<%i,%i>, (%i,%i) *= (%i,%i)", i,
            j, out.numRows(), out.numCols(), in.numRows(), in.numCols());
#endif
      elementwise_multiply<value_type, i, j>(fac, out.values(), in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(Core::LinAlg::SerialDenseMatrix::Base out,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or left.numRows() != i or left.numCols() != j or
          right.numRows() != i or right.numCols() != j)
        FOUR_C_THROW(
            "Invalid matrix sizes in elementwise_multiply<%i,%i>, (%i,%i) = (%i,%i)*(%i,%i)", i, j,
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      elementwise_multiply<value_type, i, j>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(Core::LinAlg::SerialDenseMatrix::Base out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or left.numRows() != i or left.numCols() != j or
          right.numRows() != i or right.numCols() != j)
        FOUR_C_THROW(
            "Invalid matrix sizes in elementwise_multiply<%i,%i>, (%i,%i) = (%i,%i)*(%i,%i)", i, j,
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      elementwise_multiply<value_type, i, j>(out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_multiply(const value_type outfac,
        Core::LinAlg::SerialDenseMatrix::Base out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or left.numRows() != i or left.numCols() != j or
          right.numRows() != i or right.numCols() != j)
        FOUR_C_THROW(
            "Invalid matrix sizes in elementwise_multiply<%i,%i>, (%i,%i) = (%i,%i)*(%i,%i)", i, j,
            out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      elementwise_multiply<value_type, i, j>(
          outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(
        Core::LinAlg::SerialDenseMatrix::Base out, const Core::LinAlg::SerialDenseMatrix::Base in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or in.numRows() != i or in.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_divide<%i,%i>, (%i,%i) *= (%i,%i)", i, j,
            out.numRows(), out.numCols(), in.numRows(), in.numCols());
#endif
      elementwise_divide<value_type, i, j>(out.values(), in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type fac, Core::LinAlg::SerialDenseMatrix::Base out,
        const Core::LinAlg::SerialDenseMatrix::Base in)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or in.numRows() != i or in.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_divide<%i,%i>, (%i,%i) *= (%i,%i)", i, j,
            out.numRows(), out.numCols(), in.numRows(), in.numCols());
#endif
      elementwise_divide<value_type, i, j>(fac, out.values(), in.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(Core::LinAlg::SerialDenseMatrix::Base out,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or left.numRows() != i or left.numCols() != j or
          right.numRows() != i or right.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_divide<%i,%i>, (%i,%i) = (%i,%i)*(%i,%i)",
            i, j, out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      elementwise_divide<value_type, i, j>(out.values(), left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(Core::LinAlg::SerialDenseMatrix::Base out,
        const value_type infac, const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or left.numRows() != i or left.numCols() != j or
          right.numRows() != i or right.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_divide<%i,%i>, (%i,%i) = (%i,%i)*(%i,%i)",
            i, j, out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      elementwise_divide<value_type, i, j>(out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void elementwise_divide(const value_type outfac,
        Core::LinAlg::SerialDenseMatrix::Base out, const value_type infac,
        const Core::LinAlg::SerialDenseMatrix::Base left,
        const Core::LinAlg::SerialDenseMatrix::Base right)
    {
#ifdef FOUR_C_DEBUG
      if (out.numRows() != i or out.numCols() != j or left.numRows() != i or left.numCols() != j or
          right.numRows() != i or right.numCols() != j)
        FOUR_C_THROW("Invalid matrix sizes in elementwise_divide<%i,%i>, (%i,%i) = (%i,%i)*(%i,%i)",
            i, j, out.numRows(), out.numCols(), left.numRows(), left.numCols(), right.numRows(),
            right.numCols());
#endif
      elementwise_divide<value_type, i, j>(
          outfac, out.values(), infac, left.values(), right.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void abs(Core::LinAlg::SerialDenseMatrix::Base& dest,
        const Core::LinAlg::SerialDenseMatrix::Base& src)
    {
#ifdef FOUR_C_DEBUG
      if (dest.numRows() != dest.numRows() or src.numCols() != src.numCols())
        FOUR_C_THROW("Invalid matrix sizes for abs, (%i,%i) = abs( (%i,%i) )", dest.numRows(),
            dest.numCols(), src.numRows(), src.numCols());
#endif
      abs<value_type, i, j>(dest.values(), src.values());
    }

    template <class value_type, unsigned int i, unsigned int j>
    inline void reciprocal(Core::LinAlg::SerialDenseMatrix::Base& dest,
        const Core::LinAlg::SerialDenseMatrix::Base& src)
    {
#ifdef FOUR_C_DEBUG
      if (dest.numRows() != dest.numRows() or src.numCols() != src.numCols())
        FOUR_C_THROW("Invalid matrix sizes for reciprocal, (%i,%i) = reciprocal( (%i,%i) )",
            dest.numRows(), dest.numCols(), src.numRows(), src.numCols());
#endif
      reciprocal<value_type, i, j>(dest.values(), src.values());
    }
  }  // namespace DenseFunctions


  /// Serial dense matrix with templated dimensions
  /*!
    A serial dense matrix with templated dimensions that is supposed to
    be fast and lightweight. The default scalar type is double.
    The value_type-array is allocated on the stack (small sizes, up to 512
    bytes or on the heap (larger sizes) and stored in
    column-major order, just like in Core::LinAlg::SerialDenseMatrix::Base.

    The interface is based on that of Core::LinAlg::SerialDenseMatrix::Base and
    Epetra_MultiVector. The whole View/Copy thing works a little
    different, though. See the appropriate functions for details.

    There is no operator[]. It behaves differently in
    Core::LinAlg::SerialDenseMatrix::Base and Core::LinAlg::SerialDenseVector::Base, and is not
    needed in either of them.
   */
  template <unsigned int rows, unsigned int cols, class value_type = double>
  class Matrix
  {
    static_assert(rows > 0, "Number of rows must be greater than zero!");
    static_assert(cols > 0, "Number of columns must be greater than zero!");

   private:
    /// threshold for when to allocate the memory instead of placing the matrix on the stack.
    /// set to 512 bytes (or 64 entries for double matrices).
    static constexpr bool allocatesmemory_ = rows * cols * sizeof(value_type) > 512;

    /// the pointer holding the data
    value_type* data_;

    /// for small sizes of the matrix, avoid expensive memory allocation by storing
    /// the matrix on the stack
    value_type datafieldsmall_[allocatesmemory_ ? 1 : rows * cols];

    /// whether we are a view to some other matrix
    bool isview_;

    /// only in combination with isview_. Pure read access to the underlying data.
    bool isreadonly_;

   public:
    typedef value_type scalar_type;

    /// Default constructor
    /*!
      Constructs a new Matrix and allocates the
      memory. If \e setzero==true it is filled with zeros, otherwise it
      is left uninitialized.

      \param setzero
        whether matrix should be initialised to zero
     */
    explicit Matrix(bool setzero = true);

    /// Constructor
    /*!
      Constructs a new Matrix from data \e d. If
      \e view==false (the default) the data is copied, otherwise a view to
      it is constructed.

      \param d
        pointer to data
      \param view
        whether the data is to be viewed or copied
     */
    explicit Matrix(value_type* d, bool view = false);

    /// Constructor
    /*!
      Constructs a new Matrix from data \e d. If
      \e view==false (the default) the data is copied, otherwise a view to
      it is constructed.
      \note a view is currently not possible, the data will be copied!!! a.ger 17.11.2008

      \param d
        pointer to data
      \param view
        whether the data is to be viewed or copied
     */
    explicit Matrix(const value_type* d, bool view = false);

    /// Constructor
    /*!
      Constructs a new Matrix from data \e d. If
      \e view==false (the default) the data is copied, otherwise a view to
      it is constructed.

      \param d
        matrix to be copied or viewed
      \param view
        whether the data is to be viewed or copied
     */
    explicit Matrix(Core::LinAlg::SerialDenseMatrix::Base& d, bool view = false);

    /// Constructor
    /*!
      Constructs a new Matrix from data \e d. The data is copied.

      \param d
        matrix to be copied
     */
    explicit Matrix(const Core::LinAlg::SerialDenseMatrix::Base& d);

    /// Constructor
    /*!
      Constructs a new Matrix from \e source. If
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
    Matrix(Matrix<rows, cols, value_type>& source, bool view);

    /// Copy constructor
    /*!
      Constructs a new Matrix from source. Unlike
      the Core::LinAlg::SerialDenseMatrix::Base copy constructor this one *always*
      copies the data, even when \e source is a view.

      \param source
        matrix to copy
     */
    Matrix(const Matrix<rows, cols, value_type>& source);

    /// Copy constructor
    /*!
      Constructs a new Matrix from \e source. If
      \e view==false the data is copied, otherwise a read-only view to
      it is constructed.

      \param source
        matrix to take data from
      \param view
        whether the data is to be viewed or copied

      \note This constructor sets the readonly_ flag, if \e view==true.
            In this case I recommend to use the copy constructor in combination
            with a const qualifier!
     */
    Matrix(const Matrix<rows, cols, value_type>& source, bool view);

    /// Deconstructor
    ~Matrix();

    /// Return the value_type* holding the data.
    inline const value_type* data() const { return data_; }
    /// Return the value_type* holding the data.
    inline const value_type* values() const { return data_; }
    /// Return the value_type* holding the data.
    inline value_type* data()
    {
      FOUR_C_ASSERT((not isreadonly_), "No write access to read-only data!");
      return data_;
    }
    /// Return the value_type* holding the data.
    inline value_type* values()
    {
      FOUR_C_ASSERT((not isreadonly_), "No write access to read-only data!");
      return data_;
    }
    /// Return the number of rows
    static constexpr unsigned int m() { return numRows(); }
    /// Return the number of columns
    static constexpr unsigned int n() { return numCols(); }
    /// Return the number of rows
    static constexpr unsigned int numRows() { return rows; }
    /// Return the number of columns
    static constexpr unsigned int numCols() { return cols; }
    /// Check whether the matrix is initialized
    /*!
      You cannot test whether the matrix is empty using m() and n(),
      for they will always return the templated size. Instead this
      function can be used, it tests whether the data pointer is not
      nullptr.

      \note To actually get a matrix for which is_initialized() returns
      false you must construct a view to nullptr, because the default
      constructor already allocates memory.
     */
    inline bool is_initialized() const { return data() != nullptr; }

    /// Set view
    /*!
      Set this matrix to be a view to \e data.

      \param data
        memory to be viewed
     */
    void set_view(value_type* data);

    /// Set view
    /*!
      Set this matrix to be a view to \e source.

      \param source
        matrix to be viewed
     */
    void set_view(Matrix<rows, cols, value_type>& source);

    /// Set copy
    /*!
      Set this matrix to be a copy of \e data. The difference to update(\e data)
      is that this funcion will allocate it's own memory when it was a
      view before, Update would copy the data into the view.

      \param data
        memory to copy
     */
    void set_copy(const value_type* data);

    /// Set copy
    /*!
      Set this matrix to be a copy of source. Only the value_type array
      will be copied, the \e isview_ flag is ignored (this is equivalent to
      set_copy(source.values()). The difference to update(\e source) is that this funcion will
      allocate it's own memory when it was a view before, Update would copy the data into the view.

      \param source
        matrix to copy from
     */
    void set_copy(const Matrix<rows, cols, value_type>& source);


    /// Calculate determinant
    /*!
      \return determinant
     */
    inline value_type determinant() const;

    /// invert in place
    /*!
      invert this matrix in place.

      \return determinant of matrix before inversion
     */
    inline value_type invert();

    /// invert matrix
    /*!
      invert matrix \e other and store the result in \e this.

      \param other
        matrix to be inverted
      \return determinant of \e other
     */
    inline value_type invert(const Matrix<rows, cols, value_type>& other);

    /// Set to zero
    /*!
      Sets every value in this matrix to zero. This is equivalent to
      PutScalar(0.0), but it should be faster.
     */
    inline void clear() { DenseFunctions::clear_matrix<value_type, rows, cols>(data()); }

    /// Fill with scalar
    /*!
      Sets every value in this matrix to \e scalar.

      \param scalar
        value to fill matrix with
     */
    inline void put_scalar(const value_type scalar)
    {
      DenseFunctions::put_scalar<value_type, rows, cols>(scalar, data());
    }

    /// Dot product
    /*!
      Return the dot product of \e this and \e other.

      \param other
        second factor
      \return dot product
     */
    template <class value_type_out = value_type, class value_type_other>
    inline value_type_out dot(const Matrix<rows, cols, value_type_other>& other) const
    {
      return DenseFunctions::dot<value_type_out, rows, cols>(data(), other.values());
    }

    /// Cross product
    /*!
      Return the cross product of \e left and \e right.

      \param left
      \param right
      \return cross product
     */
    template <class value_type_left, class value_type_right>
    inline void cross_product(const Matrix<rows, cols, value_type_left>& left,
        const Matrix<rows, cols, value_type_right>& right)
    {
      DenseFunctions::crossproduct<value_type, rows, cols>(data(), left.values(), right.values());
    }

    /// Compute absolute value
    /*!
      Fill this matrix with the absolute value of the numbers in \e other.

      \param other
        matrix to read values from
     */
    inline void abs(const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::abs<value_type, rows, cols>(data(), other.values());
    }

    /// Compute reciprocal value
    /*!
      Fill this matrix with the reciprocal value of the numbers in \e other.

      \param other
        matrix to read values from
     */
    inline void reciprocal(const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::reciprocal<value_type, rows, cols>(data(), other.values());
    }

    /// Scale
    /*!
      Scale matrix with \e scalar.

      \param scalar
        scaling factor
     */
    inline void scale(const value_type scalar)
    {
      DenseFunctions::scale_matrix<value_type, rows, cols>(scalar, data());
    }

    /// Copy: \e this = \e other
    /*!
      Copy \e other to \e this.

      \param other
        matrix to copy
     */
    template <class value_type_other>
    inline void update(const Matrix<rows, cols, value_type_other>& other)
    {
      DenseFunctions::update<value_type, rows, cols>(data(), other.values());
    }

    /// Scaled copy: \e this = \e scalarOther * \e other
    /*!
      Copy \e scalarOther * \e other to \e this.

      \param scalarOther
        scaling factor for other
      \param other
        matrix to read from
     */
    inline void update(const value_type scalarOther, const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::update<value_type, rows, cols>(data(), scalarOther, other.values());
    }

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
    template <class value_type_scalar_other, class value_type_other, class value_type_scalar_this>
    inline void update(const value_type_scalar_other scalarOther,
        const Matrix<rows, cols, value_type_other>& other, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::update<value_type, rows, cols>(
          scalarThis, data(), scalarOther, other.values());
    }

    /// Add: \e this = \e left + \e right
    /*!
      Store \e left + \e right in this matrix.

      \param left
        first matrix to add
      \param right
        second matrix to add
     */
    template <class value_type_left, class value_type_right>
    inline void update(const Matrix<rows, cols, value_type_left>& left,
        const Matrix<rows, cols, value_type_right>& right)
    {
      DenseFunctions::update<value_type, rows, cols>(data(), left.values(), right.values());
    }

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
    template <class value_type_scalar_left, class value_type_left, class value_type_scalar_right,
        class value_type_right>
    inline void update(const value_type_scalar_left scalarLeft,
        const Matrix<rows, cols, value_type_left>& left, const value_type_scalar_right scalarRight,
        const Matrix<rows, cols, value_type_right>& right)
    {
      DenseFunctions::update<value_type, rows, cols>(
          data(), scalarLeft, left.values(), scalarRight, right.values());
    }

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
    template <class value_type_scalar_left, class value_type_left, class value_type_scalar_right,
        class value_type_right, class value_type_scalar_this>
    inline void update(const value_type_scalar_left scalarLeft,
        const Matrix<rows, cols, value_type_left>& left, const value_type_scalar_right scalarRight,
        const Matrix<rows, cols, value_type_right>& right, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::update<value_type, rows, cols>(
          scalarThis, data(), scalarLeft, left.values(), scalarRight, right.values());
    }

    /// Transposed copy: \e this = \e other^T
    /*!
      Copy transposed \e other to \e this.

      \param other
        matrix to copy
     */
    template <class value_type_other>
    inline void update_t(const Matrix<cols, rows, value_type_other>& other)
    {
      DenseFunctions::update_t<value_type, rows, cols>(data(), other.values());
    }

    /// Scaled transposed copy: \e this = \e scalarOther * \e other^T
    /*!
      Transposed copy \e scalarOther * \e other^T to \e this.

      \param scalarOther
        scaling factor for other
      \param other
        matrix to read from
     */
    template <class value_type_other_scalar, class value_type_other>
    inline void update_t(const value_type_other_scalar scalarOther,
        const Matrix<cols, rows, value_type_other>& other)
    {
      DenseFunctions::update_t<value_type, rows, cols>(data(), scalarOther, other.values());
    }

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
    template <class value_type_other_scalar, class value_type_other, class value_type_this_scalar>
    inline void update_t(const value_type_other_scalar scalarOther,
        const Matrix<cols, rows, value_type_other>& other, const value_type_this_scalar scalarThis)
    {
      DenseFunctions::update_t<value_type, rows, cols>(
          scalarThis, data(), scalarOther, other.values());
    }

    /// Multiply element-wise: \e this(m,n) *= \e other(m,n)
    /*!
      Multiply \e this and \e other, storing the result in \e this.

      \param other
        factor
     */
    inline void elementwise_multiply(const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::elementwise_multiply<value_type, rows, cols>(data(), other.values());
    }

    /// Multiply element-wise: \e this(m,n) = \e scalar * \e this(m,n)*\e other(m,n)
    /*!
      Multiply \e this and \e other, scale by \e scalar and store the result in \e this.

      \param scalar
        scaling factor for the product
      \param other
        factor
     */
    inline void elementwise_multiply(
        const value_type scalar, const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::elementwise_multiply<value_type, rows, cols>(scalar, data(), other.values());
    }

    /// Multiply element-wise: \e this(m,n) = \e left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right and store the result in \e this.

      \param left
        first factor
      \param right
        second factor
     */
    inline void elementwise_multiply(
        const Matrix<rows, cols, value_type>& left, const Matrix<rows, cols, value_type>& right)
    {
      DenseFunctions::elementwise_multiply<value_type, rows, cols>(
          data(), left.values(), right.values());
    }

    /// Multiply element-wise: \e this(m,n) = \e scalarOther*\e left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right, scale by \e scalarOther and store the
      result in \e this.

      \param scalarOther
        scaling factor
      \param left
        first factor
      \param right
        second factor
     */
    inline void elementwise_multiply(const value_type scalarOther,
        const Matrix<rows, cols, value_type>& left, const Matrix<rows, cols, value_type>& right)
    {
      DenseFunctions::elementwise_multiply<value_type, rows, cols>(
          data(), scalarOther, left.values(), right.values());
    }

    /// Multiply element-wise: \e this(m,n) = \e scalarThis*\e this(m,n) + \e scalarOther*\e
    /// left(m,n)*\e right(m,n)
    /*!
      Multiply \e left and \e right, scale by \e scalarOther and add the result to \e this, scaled
      by \e scalarThis.

      \param scalarOther
        scaling factor the product
      \param left
        first factor, size (\c i)x(\c j)
      \param right
        second factor, size (\c i)x(\c j)
      \param scalarThis
        scaling factor for \e this
     */
    inline void elementwise_multiply(const value_type scalarOther,
        const Matrix<rows, cols, value_type>& left, const Matrix<rows, cols, value_type>& right,
        const value_type scalarThis)
    {
      DenseFunctions::elementwise_multiply<value_type, rows, cols>(
          scalarThis, data(), scalarOther, left.values(), right.values());
    }

    /// Divide element-wise: \e this(m,n) *= \e other(m,n)
    /*!
      Divide \e this by \e other, storing the result in \e this.

      \param other
        factor
     */
    inline void elementwise_divide(const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::elementwise_divide<value_type, rows, cols>(data(), other.values());
    }

    /// Divide element-wise: \e this(m,n) = \e scalar * \e this(m,n)*\e other(m,n)
    /*!
      Divide \e this by \e other, scale by \e scalar and store the result in \e this.

      \param scalar
        scaling factor for the product
      \param other
        factor
     */
    inline void elementwise_divide(
        const value_type scalar, const Matrix<rows, cols, value_type>& other)
    {
      DenseFunctions::elementwise_divide<value_type, rows, cols>(scalar, data(), other.values());
    }

    /// Divide element-wise: \e this(m,n) = \e left(m,n)*\e right(m,n)
    /*!
      Divide \e left by \e right and store the result in \e this.

      \param left
        dividend
      \param right
        divisor
     */
    inline void elementwise_divide(
        const Matrix<rows, cols, value_type>& left, const Matrix<rows, cols, value_type>& right)
    {
      DenseFunctions::elementwise_divide<value_type, rows, cols>(
          data(), left.values(), right.values());
    }

    /// Divide element-wise: \e this(m,n) = \e scalarOther*\e left(m,n)*\e right(m,n)
    /*!
      Divide \e left by \e right, scale by \e scalarOther and store the
      result in \e this.

      \param scalarOther
        scaling factor
      \param left
        dividend
      \param right
        divisor
     */
    inline void elementwise_divide(const value_type scalarOther,
        const Matrix<rows, cols, value_type>& left, const Matrix<rows, cols, value_type>& right)
    {
      DenseFunctions::elementwise_divide<value_type, rows, cols>(
          data(), scalarOther, left.values(), right.values());
    }

    /// Divide element-wise: \e this(m,n) = \e scalarThis*\e this(m,n) + \e scalarOther*\e
    /// left(m,n)/\e right(m,n)
    /*!
      Divide \e left by \e right, scale by \e scalarOther and add the result to \e this, scaled by
      \e scalarThis.

      \param scalarOther
        scaling factor the product
      \param left
        dividend, size (\c i)x(\c j)
      \param right
        divisor, size (\c i)x(\c j)
      \param scalarThis
        scaling factor for \e this
     */
    inline void elementwise_divide(const value_type scalarOther,
        const Matrix<rows, cols, value_type>& left, const Matrix<rows, cols, value_type>& right,
        const value_type scalarThis)
    {
      DenseFunctions::elementwise_divide<value_type, rows, cols>(
          scalarThis, data(), scalarOther, left.values(), right.values());
    }


    /// Calculate 1-norm
    /*!
      This is *not* the same as Core::LinAlg::SerialDenseMatrix::Base::NormOne.

      \return 1-norm
     */
    inline value_type norm1() const
    {
      return DenseFunctions::norm1<value_type, rows, cols>(data());
    }

    /// Calculate 2-norm (Euclidean norm)
    /*!
      \return 2-norm
     */
    inline value_type norm2() const
    {
      return DenseFunctions::norm2<value_type, rows, cols>(data());
    }

    /// Calculate inf-norm
    /*!
      This is *not* the same as Core::LinAlg::SerialDenseMatrix::Base::NormInf.

      \return inf-norm
     */
    inline value_type norm_inf() const
    {
      return DenseFunctions::norm_inf<value_type, rows, cols>(data());
    }

    /// Calculate minimum value
    /*!
      \return minimum value
     */
    inline value_type min_value() const
    {
      return DenseFunctions::min_value<value_type, rows, cols>(data());
    }

    /// Calculate maximum value
    /*!
      \return maximum value
     */
    inline value_type max_value() const
    {
      return DenseFunctions::max_value<value_type, rows, cols>(data());
    }

    /// Calculate mean value
    /*!
      \return mean value
     */
    inline value_type mean_value() const
    {
      return DenseFunctions::mean_value<value_type, rows, cols>(data());
    }

    /// Multiply: \e this = \e left*right
    /*!
      This is equivalent to multiply_nn(\e left,\e right).

      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_left, class value_type_right>
    inline void multiply(const Matrix<rows, inner, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right)
    {
      DenseFunctions::multiply<value_type, rows, inner, cols>(
          data(), left.values(), right.values());
    }

    /// Multiply: \e this = \e left*right
    /*!
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_left, class value_type_right>
    inline void multiply_nn(const Matrix<rows, inner, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right)
    {
      DenseFunctions::multiply<value_type, rows, inner, cols>(
          data(), left.values(), right.values());
    }

    /// Multiply: \e this = \e left*right^T
    /*!
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_left, class value_type_right>
    inline void multiply_nt(const Matrix<rows, inner, value_type_left>& left,
        const Matrix<cols, inner, value_type_right>& right)
    {
      DenseFunctions::multiply_nt<value_type, rows, inner, cols>(
          data(), left.values(), right.values());
    }

    /// Multiply: \e this = \e left^T*right
    /*!
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_left, class value_type_right>
    inline void multiply_tn(const Matrix<inner, rows, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right)
    {
      DenseFunctions::multiply_tn<value_type, rows, inner, cols>(
          data(), left.values(), right.values());
    }

    /// Multiply: \e this = \e left^T*right^T
    /*!
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_left, class value_type_right>
    inline void multiply_tt(const Matrix<inner, rows, value_type_left>& left,
        const Matrix<cols, inner, value_type_right>& right)
    {
      DenseFunctions::multiply_tt<value_type, rows, inner, cols>(
          data(), left.values(), right.values());
    }


    /// Multiply: \e this = \e scalarOthers * \e left*right
    /*!
      \param scalarOthers
        scalar factor for \e left*right
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right>
    inline void multiply(const value_type_scalar_other scalarOthers,
        const Matrix<rows, inner, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right)
    {
      DenseFunctions::multiply<value_type, rows, inner, cols>(
          data(), scalarOthers, left.values(), right.values());
    }

    /// Multiply: \e this = \e scalarOthers * \e left*right
    /*!
      This is equivalent to multiply_nn(\e scalarOthers,\e left,\e right).

      \param scalarOthers
        scalar factor for \e left*right
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right>
    inline void multiply_nn(const value_type_scalar_other scalarOthers,
        const Matrix<rows, inner, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right)
    {
      DenseFunctions::multiply<value_type, rows, inner, cols>(
          data(), scalarOthers, left.values(), right.values());
    }

    /// Multiply: \e this = \e scalarOthers * \e left*right^T
    /*!
      \param scalarOthers
        scalar factor for \e left*right^T
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right>
    inline void multiply_nt(const value_type_scalar_other scalarOthers,
        const Matrix<rows, inner, value_type_left>& left,
        const Matrix<cols, inner, value_type_right>& right)
    {
      DenseFunctions::multiply_nt<value_type, rows, inner, cols>(
          data(), scalarOthers, left.values(), right.values());
    }

    /// Multiply: \e this = \e scalarOthers * \e left^T*right
    /*!
      \param scalarOthers
        scalar factor for \e left^T*right
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right>
    inline void multiply_tn(const value_type_scalar_other scalarOthers,
        const Matrix<inner, rows, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right)
    {
      DenseFunctions::multiply_tn<value_type, rows, inner, cols>(
          data(), scalarOthers, left.values(), right.values());
    }

    /// Multiply: \e this = \e scalarOthers * \e left^T*right^T
    /*!
      \param scalarOthers
        scalar factor for \e left^T*right^T
      \param left
        first factor
      \param right
        second factor
     */
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right>
    inline void multiply_tt(const value_type_scalar_other scalarOthers,
        const Matrix<inner, rows, value_type_left>& left,
        const Matrix<cols, inner, value_type_right>& right)
    {
      DenseFunctions::multiply_tt<value_type, rows, inner, cols>(
          data(), scalarOthers, left.values(), right.values());
    }

    /// Multiply: \e this = \e scalarThis * \e this + \e scalarOthers * \e left*right
    /*!
      This is equivalent to multiply_nn(\e scalarOthers,\e left,\e right,\e scalarThis).

      \param scalarOthers
        scalar factor for \e left*right
      \param left
        first factor
      \param right
        second factor
      \param scalarThis
        scalar factor for \e this
     */
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right, class value_type_scalar_this>
    inline void multiply(const value_type_scalar_other scalarOthers,
        const Matrix<rows, inner, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::multiply<value_type, rows, inner, cols>(
          scalarThis, data(), scalarOthers, left.values(), right.values());
    }

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
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right, class value_type_scalar_this>
    inline void multiply_nn(const value_type_scalar_other scalarOthers,
        const Matrix<rows, inner, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::multiply<value_type, rows, inner, cols>(
          scalarThis, data(), scalarOthers, left.values(), right.values());
    }

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
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right, class value_type_scalar_this>
    inline void multiply_nt(const value_type_scalar_other scalarOthers,
        const Matrix<rows, inner, value_type_left>& left,
        const Matrix<cols, inner, value_type_right>& right, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::multiply_nt<value_type, rows, inner, cols>(
          scalarThis, data(), scalarOthers, left.values(), right.values());
    }

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
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right, class value_type_scalar_this>
    inline void multiply_tn(const value_type_scalar_other scalarOthers,
        const Matrix<inner, rows, value_type_left>& left,
        const Matrix<inner, cols, value_type_right>& right, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::multiply_tn<value_type, rows, inner, cols>(
          scalarThis, data(), scalarOthers, left.values(), right.values());
    }

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
    template <unsigned int inner, class value_type_scalar_other, class value_type_left,
        class value_type_right, class value_type_scalar_this>
    inline void multiply_tt(const value_type_scalar_other scalarOthers,
        const Matrix<inner, rows, value_type_left>& left,
        const Matrix<cols, inner, value_type_right>& right, const value_type_scalar_this scalarThis)
    {
      DenseFunctions::multiply_tt<value_type, rows, inner, cols>(
          scalarThis, data(), scalarOthers, left.values(), right.values());
    }

    /// Write \e this to \e out
    /*!
      Write a readable representation of \e this to \e out. This function is called by
      \e out << *\e this.

      \param out
        out stream
     */
    void print(std::ostream& out) const;

    /// = operator
    /*!
      Copy data from \e other to \e this, equivalent to update(other).
      \param other
        matrix to get data from
     */
    inline Matrix<rows, cols, value_type>& operator=(const Matrix<rows, cols, value_type>& other);

    /// = operator for double
    /*!
      Fill with double \e other, same as PutScalar(other).

      \param other
        scalar value
     */
    inline Matrix<rows, cols, value_type>& operator=(const value_type other);

    /// == operator
    /*!
      Compare \e this with \e other.

      \param other
        matrix to compare with
     */
    inline bool operator==(const Matrix<rows, cols, value_type>& other) const;

    /// != operator
    /*!
      Compare \e this with \e other.

      \param other
        matrix to compare with
     */
    inline bool operator!=(const Matrix<rows, cols, value_type>& other) const;

    /// += operator
    /*!
      Add \e other to \e this.

      \param other
        matrix to add
     */
    template <class value_type_other>
    inline Matrix<rows, cols, value_type>& operator+=(
        const Matrix<rows, cols, value_type_other>& other)
    {
      DenseFunctions::update<value_type, rows, cols>(1.0, data(), 1.0, other.values());
      return *this;
    }

    /// -= operator
    /*!
      Subtract \e other from \e this.

      \param other
        matrix to subtract
     */
    template <class value_type_other>
    inline Matrix<rows, cols, value_type>& operator-=(
        const Matrix<rows, cols, value_type_other>& other)
    {
      DenseFunctions::update<value_type, rows, cols>(1.0, data(), -1.0, other.values());
      return *this;
    }

    /// Access data
    /*!
      Return value in row \e r and column \e c.

      \param r
        row index
      \param c
        column index
     */
    inline value_type& operator()(unsigned int r, unsigned int c);

    /// Access data
    /*!
      Return value in row \e r and column \e c.

      \param r
        row index
      \param c
        column index
     */
    inline const value_type& operator()(unsigned int r, unsigned int c) const;

    /// Access data
    /*!
      Return value in row \e r. This works only for Matrices with cols==1 or rows==1 (vectors),
      otherwise a compile time error is raised.

      \param r
        index
     */
    inline value_type& operator()(unsigned int r);  // for vectors, with check at compile-time

    /// Access data
    /*!
      Return value in row \e r. This works only for Matrices with cols==1 or rows==1 (vectors),
      otherwise a compile time error is raised.

      \param r
        index
     */
    inline const value_type& operator()(unsigned int r) const;
  };

  template <class value_type, unsigned int cols, unsigned int rows>
  std::ostream& operator<<(std::ostream& out, const Matrix<rows, cols, value_type>& matrix);

  // Constructors

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(bool setzero)
      : data_(nullptr), isview_(false), isreadonly_(false)
  {
    if (allocatesmemory_)
      data_ = new value_type[rows * cols];
    else
      data_ = datafieldsmall_;
    if (setzero) DenseFunctions::clear_matrix<value_type, rows, cols>(data_);
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(value_type* d, bool view)
      : data_(nullptr), isview_(view), isreadonly_(false)
  {
    if (isview_)
    {
      data_ = d;
    }
    else
    {
      if (allocatesmemory_)
        data_ = new value_type[rows * cols];
      else
        data_ = datafieldsmall_;
      std::copy(d, d + rows * cols, data_);
    }
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(const value_type* d, bool view)
      : data_(nullptr), isview_(view), isreadonly_(false)
  {
    if (isview_)
    {
      isreadonly_ = true;
      data_ = const_cast<value_type*>(d);
    }
    else
    {
      if (allocatesmemory_)
        data_ = new value_type[rows * cols];
      else
        data_ = datafieldsmall_;
      std::copy(d, d + rows * cols, data_);
    }
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(Core::LinAlg::SerialDenseMatrix::Base& d, bool view)
      : data_(nullptr), isview_(view), isreadonly_(false)
  {
    if (d.values() == nullptr) return;
    if (d.numRows() != rows or d.numCols() != cols)
      FOUR_C_THROW("illegal matrix dimension (%d,%d)", d.numRows(), d.numCols());
    if (isview_)
    {
      data_ = d.values();
    }
    else
    {
      if (allocatesmemory_)
        data_ = new value_type[rows * cols];
      else
        data_ = datafieldsmall_;
      std::copy(d.values(), d.values() + rows * cols, data_);
    }
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(const Core::LinAlg::SerialDenseMatrix::Base& d)
      : Matrix(const_cast<Core::LinAlg::SerialDenseMatrix::Base&>(d), false)
  {
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(Matrix<rows, cols, value_type>& source, bool view)
      : data_(nullptr), isview_(view), isreadonly_(false)
  {
    if (isview_)
    {
      data_ = source.data_;
    }
    else
    {
      if (allocatesmemory_)
        data_ = new value_type[rows * cols];
      else
        data_ = datafieldsmall_;
      std::copy(source.data_, source.data_ + rows * cols, data_);
    }
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(const Matrix<rows, cols, value_type>& source)
      : data_(nullptr), isview_(false), isreadonly_(false)
  {
    if (allocatesmemory_)
      data_ = new value_type[rows * cols];
    else
      data_ = datafieldsmall_;
    std::copy(source.data_, source.data_ + rows * cols, data_);
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::Matrix(const Matrix<rows, cols, value_type>& source, bool view)
      : data_(nullptr), isview_(view), isreadonly_(false)
  {
    if (isview_)
    {
      isreadonly_ = true;
      data_ = const_cast<value_type*>(source.values());
    }
    else
    {
      if (allocatesmemory_)
        data_ = new value_type[rows * cols];
      else
        data_ = datafieldsmall_;
      std::copy(source.data_, source.data_ + rows * cols, data_);
    }
  }

  // Destructor
  template <unsigned int rows, unsigned int cols, class value_type>
  Matrix<rows, cols, value_type>::~Matrix()
  {
    if (allocatesmemory_ && not isview_) delete[] data_;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  void Matrix<rows, cols, value_type>::set_view(value_type* data)
  {
    if (not isview_)
    {
      if (allocatesmemory_) delete[] data_;
      isview_ = true;
    }
    data_ = data;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  void Matrix<rows, cols, value_type>::set_view(Matrix<rows, cols, value_type>& source)
  {
    if (not isview_)
    {
      if (allocatesmemory_) delete[] data_;
      isview_ = true;
    }
    data_ = source.data_;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  void Matrix<rows, cols, value_type>::set_copy(const value_type* data)
  {
    if (isview_)
    {
      if (allocatesmemory_)
        data_ = new value_type[rows * cols];
      else
        data_ = datafieldsmall_;
      isview_ = false;
    }
    std::copy(data, data + rows * cols, data_);
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  void Matrix<rows, cols, value_type>::set_copy(const Matrix<rows, cols, value_type>& source)
  {
    set_copy(source.values());
  }

  // Determinant
  template <unsigned int rows, unsigned int cols, class value_type>
  inline value_type Matrix<rows, cols, value_type>::determinant() const
  {
    static_assert(rows == cols, "Cannot compute determinant of non-square matrix.");
    return DenseFunctions::determinant<value_type, rows, cols>(data());
  }

  // invert
  template <unsigned int rows, unsigned int cols, class value_type>
  inline value_type Matrix<rows, cols, value_type>::invert()
  {
    static_assert(rows == cols, "Cannot compute inverse of non-square matrix");
    return DenseFunctions::invert<value_type, rows, cols>(data());
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline value_type Matrix<rows, cols, value_type>::invert(
      const Matrix<rows, cols, value_type>& other)
  {
    static_assert(rows == cols, "Cannot compute inverse of non-square matrix");
    return DenseFunctions::invert<value_type, rows, cols>(data(), other.values());
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  void Matrix<rows, cols, value_type>::print(std::ostream& out) const
  {
    out << "Matrix<" << rows << ',' << cols << '>';
    if (isview_) out << " (view to memory only)";
    if (isreadonly_) out << "\n (read only)";
    if (data() == nullptr)
    {
      out << " with data_==nullptr!\n";
      return;
    }
    if (cols > 1)
    {
      out << "\n[";
      for (unsigned int i = 0; i < rows; ++i)
      {
        if (i != 0) out << ' ';
        for (unsigned int j = 0; j < cols; ++j)
        {
          out << data()[i + rows * j];
          if (j + 1 < cols) out << ", ";
        }
        if (i + 1 < rows)
          out << ",\n";
        else
          out << "]\n";
      }
    }
    else
    {
      out << "[";
      for (unsigned int i = 0; i < rows; ++i)
      {
        if (i != 0) out << ' ';
        out << data()[i];
      }
      out << "]\n";
    }
  }

  /// output operator for Matrix
  /*!
    Write matrix to out. This function calls matrix.print(out).
   */
  template <class value_type, unsigned int cols, unsigned int rows>
  std::ostream& operator<<(std::ostream& out, const Matrix<rows, cols, value_type>& matrix)
  {
    matrix.print(out);
    return out;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline Matrix<rows, cols, value_type>& Matrix<rows, cols, value_type>::operator=(
      const Matrix<rows, cols, value_type>& other)
  {
    DenseFunctions::update<value_type, rows, cols>(data(), other.values());
    return *this;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline Matrix<rows, cols, value_type>& Matrix<rows, cols, value_type>::operator=(
      const value_type other)
  {
    DenseFunctions::put_scalar<value_type, rows, cols>(other, data());
    return *this;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline bool Matrix<rows, cols, value_type>::operator==(
      const Matrix<rows, cols, value_type>& other) const
  {
    if (data() == other.values()) return true;
    // unfortunately memcmp does not work, because +0 and -0 are
    // different in memory...
    for (unsigned int c = 0; c < rows * cols; ++c)
      if (data()[c] != other.values()[c]) return false;
    return true;
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline bool Matrix<rows, cols, value_type>::operator!=(
      const Matrix<rows, cols, value_type>& other) const
  {
    return not(*this == other);
  }

  // Access operator

  template <unsigned int rows, unsigned int cols, class value_type>
  inline value_type& Matrix<rows, cols, value_type>::operator()(unsigned int r, unsigned int c)
  {
#ifdef FOUR_C_DEBUG
    if (r >= rows or c >= cols)
      FOUR_C_THROW("Indices %i,%i out of range in Matrix<%i,%i>.", r, c, rows, cols);
#endif
    return data()[r + c * rows];
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline const value_type& Matrix<rows, cols, value_type>::operator()(
      unsigned int r, unsigned int c) const
  {
#ifdef FOUR_C_DEBUG
    if (r >= rows or c >= cols)
      FOUR_C_THROW("Indices %i,%i out of range in Matrix<%i,%i>.", r, c, rows, cols);
#endif
    return data()[r + c * rows];
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline value_type& Matrix<rows, cols, value_type>::operator()(unsigned int r)
  {
    static_assert((cols == 1) or (rows == 1), "cannot call 1-d access function on 2-d matrix");
#ifdef FOUR_C_DEBUG
    if (r >= (cols == 1 ? rows : cols))
      FOUR_C_THROW("Index %i out of range in Matrix<%i,%i>.", r, rows, cols);
#endif
    return data()[r];
  }

  template <unsigned int rows, unsigned int cols, class value_type>
  inline const value_type& Matrix<rows, cols, value_type>::operator()(unsigned int r) const
  {
    static_assert((cols == 1) or (rows == 1), "cannot call 1-d access function on 2-d matrix");
#ifdef FOUR_C_DEBUG
    if (r >= (cols == 1 ? rows : cols))
      FOUR_C_THROW("Index %i out of range in Matrix<%i,%i>.", r, rows, cols);
#endif
    return data()[r];
  }


  /// A solver for fixed size serial dense matrices
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
    void SetMatrix(Matrix<rows, cols, double>& mat);

    /// Set the vectors
    /*!
      Set the vectors, the new equation is matrix*X=B.

      \param X
        vector of unknowns
      \param B
        right hand side vector
     */
    void SetVectors(Matrix<cols, dim_rhs, double>& X, Matrix<rows, dim_rhs, double>& B);

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
    void SolveWithTranspose(bool b) { transpose_ = b; }

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
  void FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::SetMatrix(Matrix<rows, cols, double>& mat)
  {
    c_.clear();
    r_.clear();
    pivot_vec_.clear();
    inverted_ = factored_ = solved_ = false;
    // vec_B_ = vec_X_ = nullptr;
    matrix_ = &mat;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  void FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::SetVectors(
      Matrix<cols, dim_rhs, double>& X, Matrix<rows, dim_rhs, double>& B)
  {
    solved_ = false;
    vec_x_ = &X;
    vec_b_ = &B;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  void FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::factor_with_equilibration(bool b)
  {
#ifdef FOUR_C_DEBUG
    if (factored_ or inverted_)
      FOUR_C_THROW("Cannot set equilibration after changing the matrix with Factor() or invert().");
#endif
    equilibrate_ = b;
  }

  template <unsigned int rows, unsigned int cols, unsigned int dim_rhs>
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::Factor()
  {
#ifdef FOUR_C_DEBUG
    if (inverted_) FOUR_C_THROW("Cannot factor the inverted matrix.");
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
  int FixedSizeSerialDenseSolver<rows, cols, dim_rhs>::Solve()
  {
    int errnum = 0;
    if (equilibrate_)
    {
      errnum = equilibrate_rhs();
    }
    if (errnum != 0) return errnum;
#ifdef FOUR_C_DEBUG
    if (not vec_b_ or not vec_x_) FOUR_C_THROW("Both vectors must be set to solve.");
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
        errnum = Factor();
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
      errnum = Factor();
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
