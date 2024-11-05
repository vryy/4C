// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_SSTI_HPP
#define FOUR_C_INPAR_SSTI_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class ConditionDefinition;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace SSTI
  {
    /// Type of coupling strategy for SSI problems
    enum class SolutionScheme
    {
      monolithic
    };

    //! type of scalar transport time integration
    enum class ScaTraTimIntType
    {
      elch
    };

    /// set the ssti parameters
    void set_valid_parameters(Teuchos::ParameterList& list);

    /// set specific ssti conditions
    void set_valid_conditions(
        std::vector<std::shared_ptr<Core::Conditions::ConditionDefinition>>& condlist);

  }  // namespace SSTI
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
