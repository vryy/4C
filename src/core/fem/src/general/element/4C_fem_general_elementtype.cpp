// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_elementtype.hpp"

FOUR_C_NAMESPACE_OPEN

Core::Elements::ElementType::ElementType() : ParObjectType() {}

int Core::Elements::ElementType::initialize(Core::FE::Discretization& dis) { return 0; }

FOUR_C_NAMESPACE_CLOSE
