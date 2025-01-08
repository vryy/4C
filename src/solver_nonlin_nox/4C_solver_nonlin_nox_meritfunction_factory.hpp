// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_MERITFUNCTION_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_MERITFUNCTION_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_globaldata.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace MeritFunction
    {
      class Factory
      {
       public:
        //! constructor
        Factory();

        //! destructor
        virtual ~Factory() = default;

        /** \brief get a valid merit function pointer
         *
         *  choose between constraint and unconstraint */
        Teuchos::RCP<::NOX::MeritFunction::Generic> build_merit_function(
            const NOX::Nln::GlobalData& noxNlnGlobalData) const;

       private:
        //! unconstraint factory
        Teuchos::RCP<::NOX::MeritFunction::Generic> build_unconstrained_merit_function(
            const NOX::Nln::MeritFunction::MeritFctName& mftype,
            const NOX::Nln::GlobalData& noxNlnGlobalData) const;

        //! constraint factory
        Teuchos::RCP<::NOX::MeritFunction::Generic> build_constrained_merit_function(
            const NOX::Nln::MeritFunction::MeritFctName& mftype,
            const NOX::Nln::GlobalData& noxNlnGlobalData) const;

      };  // class Factory

      /*! \brief Non-member function to build a merit function object.

      \relates NOX::NlnSol::Constraint::MeritFunction::Factory

      */
      Teuchos::RCP<::NOX::MeritFunction::Generic> build_merit_function(
          const NOX::Nln::GlobalData& noxNlnGlobalData);
    }  // namespace MeritFunction
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
