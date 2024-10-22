// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_ART_TERMINAL_BC_HPP
#define FOUR_C_ART_NET_ART_TERMINAL_BC_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN



namespace Arteries
{
  namespace Utils
  {
    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for inlet bc
    //--------------------------------------------------------------------

    /*!
    \brief 1d-artery inlet boundary condition, this class is meant to do
     solve the bc at the inlet of a one-dimensional arterial network
    */

    /*
    \brief Standard solver for 1d artery inlet
    */
    void solve_prescribed_terminal_bc(Core::FE::Discretization& actdis,
        const Core::Conditions::Condition* condition, Teuchos::ParameterList& params);

    /*
    \brief Standard solver for 1d artery reflective outlet
    */
    void solve_reflective_terminal(Core::FE::Discretization& actdis,
        const Core::Conditions::Condition* condition, Teuchos::ParameterList& params);

    /*
    \brief Standard solver for 1d artery explicit windkessel BC outlet
    */
    void solve_expl_windkessel_bc(Core::FE::Discretization& actdis,
        const Core::Conditions::Condition* condition, Teuchos::ParameterList& params);

  }  // namespace Utils
}  // namespace Arteries

FOUR_C_NAMESPACE_CLOSE

#endif
