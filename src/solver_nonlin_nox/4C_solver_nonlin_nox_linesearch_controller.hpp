// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_CONTROLLER_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_CONTROLLER_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace LineSearch
    {
      class Controller
      {
       public:
        //! @name Access functionality
        //@{
        //! get the number of line search iterations
        [[nodiscard]] virtual int get_num_iterations() const = 0;

        //! get the merit function
        [[nodiscard]] virtual const ::NOX::MeritFunction::Generic& get_merit_function() const = 0;

        //! get the current search direction
        [[nodiscard]] virtual const ::NOX::Abstract::Vector& get_search_direction() const = 0;

        //! get current step length
        [[nodiscard]] virtual double get_step_length() const = 0;

        //!@}

        //! @name Mutator functionality
        //! @{
        //! set current step length
        virtual void set_step_length(double step) = 0;
        //! @}
      };
    }  // namespace LineSearch
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
