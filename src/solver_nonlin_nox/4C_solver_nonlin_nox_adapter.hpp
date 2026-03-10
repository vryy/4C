// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_ADAPTER_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_ADAPTER_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_interface_required_base.hpp"

#include <mpi.h>
#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Solver_Generic.H>
#include <NOX_StatusTest_Generic.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <functional>
#include <map>
#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
  class SparseOperator;
  template <typename T>
  class Vector;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }
  namespace Nln
  {
    class GlobalData;
    namespace Inner::StatusTest
    {
      class Generic;
    }
    class LinearSystemBase;
    class Problem;
    class Vector;
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest

    /**
     * @brief Convenience adapter that wraps NOX behind callbacks.
     *
     * Users provide residual/Jacobian callbacks and bound jacobians in the
     * constructor, then call @ref solve. The adapter manages NOX internals (GlobalData, Problem,
     * LinearSystem, Group, StatusTests, Solver).
     *
     * Example:
     * @code
     *   std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>> solvers;
     *   solvers[NOX::Nln::sol_generic] = linear_solver;
     *
     *   auto residual_callback = [&](const Core::LinAlg::Vector<double>& x,
     *                       Core::LinAlg::Vector<double>& f,
     *                       NOX::Nln::FillType) {
     *     // fill residual f(x)
     *     return true;
     *   };
     *   auto jacobian_callback = [&](const Core::LinAlg::Vector<double>& x,
     *                                Core::LinAlg::SparseOperator& J) {
     *     // assemble Jacobian J(x)
     *     return true;
     *   };
     *
     *   NOX::Nln::Adapter adapter(
     *       comm, nox_params, solvers, x, jacobian, residual_callback,
     * jacobian_callback);
     *
     *   const unsigned int iterations = adapter.solve();
     * @endcode
     */
    class Adapter
    {
     public:
      using ResidualCallback = std::function<bool(
          const Core::LinAlg::Vector<double>&, Core::LinAlg::Vector<double>&, NOX::Nln::FillType)>;
      using JacobianCallback =
          std::function<bool(const Core::LinAlg::Vector<double>&, Core::LinAlg::SparseOperator&)>;

      /**
       * @brief Construct adapter and build the NOX stack.
       *
       * @param comm MPI communicator.
       * @param nox_params NOX parameter list.
       * @param solvers Map of linear solvers keyed by solution type.
       * @param x In/out solution vector bound to this adapter instance.
       * @param jacobian Jacobian operator bound to this adapter instance.
       * @param residual_callback Residual callback.
       * @param jacobian_callback Jacobian callback.
       */
      Adapter(MPI_Comm comm, const Teuchos::ParameterList& nox_params,
          const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
          Core::LinAlg::Vector<double>& x, Core::LinAlg::SparseOperator& jacobian,
          ResidualCallback residual_callback, JacobianCallback jacobian_callback);

      /**
       * @brief Run nonlinear solve.
       *
       * @return Number of nonlinear iterations performed.
       */
      unsigned int solve();

     private:
      /// Create the full NOX problem/group/solver stack for the bound vectors/operators.
      void build_solver();

      Teuchos::ParameterList nox_params_;
      std::shared_ptr<NOX::Nln::Vector> x_nox_;
      std::shared_ptr<Core::LinAlg::SparseOperator> jacobian_;
      Teuchos::RCP<NOX::Nln::GlobalData> nox_global_data_;
      Teuchos::RCP<NOX::Nln::Problem> nox_problem_;
      Teuchos::RCP<NOX::Nln::LinearSystemBase> linsys_;
      Teuchos::RCP<::NOX::Abstract::Group> group_;
      Teuchos::RCP<::NOX::StatusTest::Generic> ostatus_;
      Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> istatus_;
      Teuchos::RCP<::NOX::Solver::Generic> nox_solver_;
    };
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
