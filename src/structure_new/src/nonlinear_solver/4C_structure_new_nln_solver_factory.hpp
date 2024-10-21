// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace Nln
  {
    namespace SOLVER
    {
      class Generic;
      /*! \brief Factory to build the desired non-linear structural solver
       *
       *  \author Michael Hiermeier */
      class Factory
      {
       public:
        //! constructor
        Factory();

        //! destructor
        virtual ~Factory() = default;

        //! build the specific nonlinear solver
        Teuchos::RCP<Solid::Nln::SOLVER::Generic> build_nln_solver(
            const enum Inpar::Solid::NonlinSolTech& nlnSolType) const;
      };

      /*! Non-member function, which relates to the Solid::Nln::SOLVER::Factory class
       *  Please call this method from outside! */
      Teuchos::RCP<Solid::Nln::SOLVER::Generic> build_nln_solver(
          const enum Inpar::Solid::NonlinSolTech& nlnSolType);
    }  // namespace SOLVER
  }    // namespace Nln
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
