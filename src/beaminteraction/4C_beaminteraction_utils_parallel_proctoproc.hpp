// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_UTILS_PARALLEL_PROCTOPROC_HPP
#define FOUR_C_BEAMINTERACTION_UTILS_PARALLEL_PROCTOPROC_HPP


#include "4C_config.hpp"

#include <map>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Utils
  {
    //! send data to rank map key and recv data
    void i_send_receive_any(Core::FE::Discretization& discret,
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosenddata,
        std::vector<std::pair<int, std::vector<int>>>& recvdata);

  }  // namespace Utils
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
