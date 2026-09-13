// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_SERIALDENSEVECTOR_HPP
#define FOUR_C_LINALG_SERIALDENSEVECTOR_HPP


#include "4C_config.hpp"

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <ostream>
#include <span>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * \brief Composition wrapper for Teuchos::SerialDenseVector.
   *
   * Access to the underlying Trilinos object is explicit via base().
   */
  class FOUR_C_API(FOUR_C_CORE) SerialDenseVector
  {
   public:
    using ordinalType = int;
    using scalarType = double;
    using Base = Teuchos::SerialDenseVector<ordinalType, scalarType>;

    /** \name Construction */
    //@{
    SerialDenseVector() = default;
    SerialDenseVector(int length);
    SerialDenseVector(int length, bool zeroOut);

    /*! \brief Construct from user-provided storage.
     *
     * Behavior (view/copy) follows \p cv semantics of Teuchos::DataAccess.
     */
    SerialDenseVector(Teuchos::DataAccess cv, double* values, int length);

    /*! \brief Construct from underlying Trilinos vector.
     *
     * \note \p trans is accepted for API compatibility and validated.
     */
    SerialDenseVector(const Base& source, Teuchos::ETransp trans = Teuchos::NO_TRANS);
    //@}

    // NOLINTBEGIN(readability-identifier-naming)

    /** \name Size queries and norms */
    //@{
    [[nodiscard]] int length() const;
    [[nodiscard]] int num_rows() const;
    [[nodiscard]] int numRows() const;
    [[nodiscard]] double normInf() const;
    [[nodiscard]] double normOne() const;
    [[nodiscard]] double norm2() const;
    [[nodiscard]] double normFrobenius() const;
    [[nodiscard]] bool empty() const;
    //@}

    /** \name Element access */
    //@{
    double& operator()(int i) { return vec_(i); }
    const double& operator()(int i) const { return vec_(i); }
    double& operator[](int i) { return vec_[i]; }
    const double& operator[](int i) const { return vec_[i]; }
    //@}

    /** \name Raw data access */
    //@{
    /*! \brief Returns pointer to contiguous vector storage.
     *
     * The const overload follows Teuchos semantics and returns mutable data.
     */
    double* values() const;

    //! Alias for values().
    double* data() const;
    //@}

    /** \name Modifiers */
    //@{
    int size(int length);
    int Size(int length);
    int resize(int length);
    int scale(double alpha);
    int put_scalar(double val = 0.0);
    int putScalar(double val = 0.0);
    void assign(const SerialDenseVector& source);
    //@}

    /** \name Algebraic operations */
    //@{
    SerialDenseVector& operator+=(const SerialDenseVector& other);
    double dot(const SerialDenseVector& other) const;
    int multiply(Teuchos::ETransp transa, Teuchos::ETransp transb, double alpha,
        const Teuchos::SerialDenseMatrix<ordinalType, scalarType>& A,
        const Teuchos::SerialDenseMatrix<ordinalType, scalarType>& B, double beta);
    int multiply(Teuchos::ETransp transa, Teuchos::ETransp transb, double alpha,
        const Teuchos::SerialDenseMatrix<ordinalType, scalarType>& A, const SerialDenseVector& B,
        double beta);
    void print(std::ostream& out) const;
    //@}

    // NOLINTEND(readability-identifier-naming)

    /** \name Trilinos interoperability */
    //@{
    /*! \brief Access the wrapped Trilinos vector for low-level interfaces. */
    Base& base();

    /*! \brief Const access variant of base(). */
    const Base& base() const;
    //@}

    /** \name Span views */
    //@{
    //! Returns a read-only span view over the current vector data.
    [[nodiscard]] std::span<const double> as_span() const { return std::span(values(), length()); }

    //! Returns a mutable span view over the current vector data.
    [[nodiscard]] std::span<double> as_span() { return std::span(values(), length()); }
    //@}

   private:
    /*! \cond INTERNAL */
    Base vec_;
    /*! \endcond */
  };

  // type definition for serial integer vector
  using IntSerialDenseVector = Teuchos::SerialDenseVector<int, int>;

  /** \name Free vector utilities */
  //@{

  /*!
    \brief Update vector components with scaled values.

    Computes \p b = \p alpha * \p a + \p beta * \p b.
    */
  void update(double alpha, const SerialDenseVector& a, double beta, SerialDenseVector& b);

  //! Euclidean vector norm.
  double norm2(const SerialDenseVector& v);

  //@}

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const SerialDenseVector& vec)
  {
    vec.print(out);
    return out;
  }

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const IntSerialDenseVector& vec)
  {
    vec.print(out);
    return out;
  }
}  // namespace Core::LinAlg


FOUR_C_NAMESPACE_CLOSE

#endif
