// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_IO_MONITOR_STRUCTURE_DBC_HPP
#define FOUR_C_INPAR_IO_MONITOR_STRUCTURE_DBC_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace IOMonitorStructureDBC
  {
    /// data format for written numeric data
    enum FileType
    {
      csv,
      data
    };

    /// set the valid parameters related to writing of output at runtime
    void set_valid_parameters(Teuchos::ParameterList& list);

  }  // namespace IOMonitorStructureDBC
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
