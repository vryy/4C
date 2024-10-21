// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_DG_ELEMENT_HPP
#define FOUR_C_FEM_GENERAL_DG_ELEMENT_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class DgElement
  {
   public:
    /*!
     \brief Destructor
    */
    virtual ~DgElement() = default;
    virtual int num_dof_per_node_auxiliary() const = 0;

    virtual int num_dof_per_element_auxiliary() const = 0;
  };
}  // namespace Core::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
