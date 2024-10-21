// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SOLVER_PREPOSTOP_GENERIC_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SOLVER_PREPOSTOP_GENERIC_HPP

#include "4C_config.hpp"

#include <NOX_Observer.hpp>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace Solver
    {
      namespace PrePostOp
      {
        /*! \brief Non-linear solver helper class
         *
         * The NOX::Nln::Solver s are supposed to use this helper class to save
         * solver dependent things, such as the number of nonlinear solver steps.
         *
         * ToDo Extend the functionality to non-linear solvers which are not
         * derived from the LineSearchBased class.
         *
         * \author Michael Hiermeier */

        class Generic : public ::NOX::Observer
        {
         public:
          //! constructor
          Generic();

          void runPreIterate(const ::NOX::Solver::Generic& solver) override;

          void runPreSolve(const ::NOX::Solver::Generic& nlnSolver) override;
        };
      }  // namespace PrePostOp
    }    // namespace Solver
  }      // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
