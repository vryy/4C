// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_NOX_AITKEN_HPP
#define FOUR_C_FSI_NOX_AITKEN_HPP

#include "4C_config.hpp"

#include <NOX_GlobalData.H>
#include <NOX_LineSearch_Generic.H>  // base class
#include <NOX_LineSearch_UserDefinedFactory.H>
#include <NOX_Utils.H>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace FSI
  {
    //! Aikten line search - the simple relaxation.
    /*!
      This line search can be called via ::NOX::LineSearch::Manager.

      The working horse in FSI.

     */
    class AitkenRelaxation : public ::NOX::LineSearch::Generic
    {
     public:
      //! Constructor
      AitkenRelaxation(const Teuchos::RCP<::NOX::Utils>& utils, Teuchos::ParameterList& params);


      // derived
      bool reset(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);

      // derived
      bool compute(::NOX::Abstract::Group& newgrp, double& step, const ::NOX::Abstract::Vector& dir,
          const ::NOX::Solver::Generic& s) override;

      //! return relaxation parameter
      double get_omega();

     private:
      //! difference of last two solutions
      Teuchos::RCP<::NOX::Abstract::Vector> del_;

      //! difference of difference of last two pair of solutions
      Teuchos::RCP<::NOX::Abstract::Vector> del2_;

      //! aitken factor
      double nu_;

      //! max step size
      double maxstep_;

      //! min step size
      double minstep_;

      //! flag for restart
      bool restart_;

      //! first omega after restart
      double restart_omega_;

      //! Printing utilities
      Teuchos::RCP<::NOX::Utils> utils_;
    };


    /// simple factory that creates aitken relaxation class
    class AitkenFactory : public ::NOX::LineSearch::UserDefinedFactory
    {
     public:
      Teuchos::RCP<::NOX::LineSearch::Generic> buildLineSearch(
          const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) const override
      {
        if (aitken_ == Teuchos::null)
          aitken_ = Teuchos::make_rcp<AitkenRelaxation>(gd->getUtils(), params);
        else
          aitken_->reset(gd, params);
        return aitken_;
      }

      Teuchos::RCP<AitkenRelaxation> get_aitken() { return aitken_; };

     private:
      mutable Teuchos::RCP<AitkenRelaxation> aitken_;
    };

  }  // namespace FSI
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
