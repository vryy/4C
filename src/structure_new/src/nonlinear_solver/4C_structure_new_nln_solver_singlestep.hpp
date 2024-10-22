// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_SINGLESTEP_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_SINGLESTEP_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_statustest_factory.hpp"
#include "4C_structure_new_nln_solver_nox.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid::Nln::SOLVER
{
  /*! \brief Full step with single iteration via NOX for explicit structural dynamics
   *
   */
  class SingleStep : public Nox
  {
   public:
    //! derived from the base class
    void setup() override;

    //! derived from the base class
    Inpar::Solid::ConvergenceStatus solve() override;

   protected:
    //! Reset the non-linear solver parameters and variables
    void reset_params() override;

   private:
    //! set the single step parameters in the nox parameter list
    void set_single_step_params();

    //! set the single step parameters from the parameter list
    void set_single_step_params(Teuchos::ParameterList& p);
  };  // class SingleStep
}  // namespace Solid::Nln::SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
