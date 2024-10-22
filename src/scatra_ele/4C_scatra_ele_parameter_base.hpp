// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_PARAMETER_BASE_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_BASE_HPP


#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

// forward declaration
FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterBase
    {
     public:
      /// Virtual destructor.
      virtual ~ScaTraEleParameterBase() = default;

      //! set parameters
      virtual void set_parameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
