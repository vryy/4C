// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_SERIALDENSESOLVER_HPP
#define FOUR_C_LINALG_SERIALDENSESOLVER_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * \brief Wrapper for Teuchos::SerialDenseSolver.
   *
   * \note The solver stores non-owning references to matrix/vector arguments
   *       via Teuchos::RCP. Referenced objects must remain alive until
   *       factor(), solve(), and invert() complete.
   */
  class FOUR_C_API(FOUR_C_CORE) SerialDenseSolver
  {
   public:
    using Base = Teuchos::SerialDenseSolver<int, double>;

    /** \name System setup */
    //@{
    /*! \brief Set system matrix A of the system A*x = b. */
    void set_matrix(SerialDenseMatrix& matrix);

    /*!
     * \brief Set vector unknowns and right-hand side for A*x = b.
     *
     * \param lhs Unknown vector x (solution written here).
     * \param rhs Right-hand side vector b.
     */
    void set_vectors(SerialDenseVector& lhs, SerialDenseVector& rhs);

    /*!
     * \brief Set matrix unknowns and right-hand side for multi-RHS systems of A*x = b.
     *
     * \param lhs Unknown x (stored as a matrix; solution written here).
     * \param rhs Right-hand side b (stored as a matrix).
     */
    void set_vectors(SerialDenseMatrix& lhs, SerialDenseMatrix& rhs);
    //@}

    /** \name Solver options */
    //@{
    /*! \brief Enable or disable matrix equilibration before factorization. */
    void factor_with_equilibration(bool enable);

    /*! \brief Enable or disable iterative refinement in solve(). */
    void solve_to_refined_solution(bool enable);
    //@}

    /** \name Operations */
    //@{
    /*!
     * \brief Factor current matrix.
     * \return Trilinos/LAPACK-style error code (0 on success).
     */
    int factor();

    /*!
     * \brief Solve for current vectors with current matrix settings.
     * \return Trilinos/LAPACK-style error code (0 on success).
     */
    int solve();

    /*!
     * \brief Invert current matrix in place.
     * \return Trilinos/LAPACK-style error code (0 on success).
     */
    int invert();
    //@}

   private:
    Base solver_;
  };
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
