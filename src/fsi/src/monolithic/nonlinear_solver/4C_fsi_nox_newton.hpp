// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_NOX_NEWTON_HPP
#define FOUR_C_FSI_NOX_NEWTON_HPP

#include "4C_config.hpp"

#include "4C_inpar_fsi.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <NOX_Direction_Newton.H>
#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace FSI
  {
    /// NOX Newton direction with adaptive linear solver tolerance
    class Newton : public ::NOX::Direction::Newton
    {
     public:
      Newton(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);


      // derived
      bool reset(
          const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) override;

      // derived
      bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
          const ::NOX::Solver::Generic& solver) override;

      void residual(double current, double desired);

     private:
      /// Printing Utilities
      Teuchos::RCP<::NOX::Utils> utils_;

      //! "Direction" sublist with parameters for the direction vector
      Teuchos::ParameterList* params_ptr_;

      /// nonlinear tolerance we strive to achieve
      double desirednlnres_;

      /// current nonlinear tolerance (what we gained after the last linear solve)
      double currentnlnres_;

      /// basic (unmodified) linear solver (AZ_r0) tolerance
      double plaintol_;

      /// improvement factor
      double better_;

      /// verbosity level of FSI algorithm
      Inpar::FSI::Verbosity verbosity_;

      std::vector<double> cresiduals_;
      std::vector<double> dresiduals_;
    };
  }  // namespace FSI
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
