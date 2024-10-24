// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SHELL7P_ELE_FACTORY_HPP
#define FOUR_C_SHELL7P_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_inpar_structure.hpp"

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  // forward declaration
  class Shell7pEleCalcInterface;
  class Shell7p;

  class Shell7pFactory
  {
   public:
    //! copy constructor
    Shell7pFactory() = default;

    static std::unique_ptr<Shell7pEleCalcInterface> provide_shell7p_calculation_interface(
        const Core::Elements::Element& ele, const std::set<Inpar::Solid::EleTech>& eletech);

   private:
    //! define shell calculation instances dependent on element technology
    template <Core::FE::CellType distype>
    static std::unique_ptr<Shell7pEleCalcInterface> define_calculation_interface_type(
        const std::set<Inpar::Solid::EleTech>& eletech);
  };  // class Shell7pFactory
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
