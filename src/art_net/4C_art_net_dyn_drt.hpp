// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_DYN_DRT_HPP
#define FOUR_C_ART_NET_DYN_DRT_HPP


#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Adapter
{
  class ArtNet;
}

void dyn_art_net_drt();

Teuchos::RCP<Adapter::ArtNet> dyn_art_net_drt(bool CoupledTo3D);


FOUR_C_NAMESPACE_CLOSE

#endif
