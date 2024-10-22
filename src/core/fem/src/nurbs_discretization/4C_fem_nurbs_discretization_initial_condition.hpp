// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_NURBS_DISCRETIZATION_INITIAL_CONDITION_HPP
#define FOUR_C_FEM_NURBS_DISCRETIZATION_INITIAL_CONDITION_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::LinAlg
{
  class Solver;
}

namespace Core::FE
{
  namespace Nurbs
  {
    /*----------------------------------------------------------------------*/
    /*!
    \brief A service method allowing the application of initial conditions
           for nurbs discretisations. Recommended version with separate
           solver allocation

    \param dis            (i) the discretisation
    \param solverparams   (i) a list with solver parameters
    \param start_function (i) a function defining the initial field (i.e. u_0(x))
    \param initialvals    (o) the initial field on output (i.e. u_cp)

    \date 08/11
    */
    void apply_nurbs_initial_condition(Core::FE::Discretization& dis,
        const Teuchos::ParameterList& solverparams,
        const Core::Utils::FunctionOfSpaceTime& start_function,
        Teuchos::RCP<Core::LinAlg::Vector<double>> initialvals);

  }  // namespace Nurbs

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
