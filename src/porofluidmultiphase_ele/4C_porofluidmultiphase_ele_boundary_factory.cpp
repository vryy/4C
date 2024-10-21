// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluidmultiphase_ele_boundary_factory.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_porofluidmultiphase_ele_boundary_calc.hpp"
#include "4C_porofluidmultiphase_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN


/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseEleInterface*
Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::provide_impl(
    const Core::Elements::Element* ele, const int numdofpernode, const std::string& disname)
{
  switch (ele->shape())
  {
    case Core::FE::CellType::quad4:
    {
      return define_problem_type<Core::FE::CellType::quad4>(numdofpernode, disname);
    }
    case Core::FE::CellType::quad8:
    {
      return define_problem_type<Core::FE::CellType::quad8>(numdofpernode, disname);
    }
    case Core::FE::CellType::quad9:
    {
      return define_problem_type<Core::FE::CellType::quad9>(numdofpernode, disname);
    }
    case Core::FE::CellType::tri3:
    {
      return define_problem_type<Core::FE::CellType::tri3>(numdofpernode, disname);
    }
    case Core::FE::CellType::tri6:
    {
      return define_problem_type<Core::FE::CellType::tri6>(numdofpernode, disname);
    }
    case Core::FE::CellType::line2:
    {
      return define_problem_type<Core::FE::CellType::line2>(numdofpernode, disname);
    }
    case Core::FE::CellType::line3:
    {
      return define_problem_type<Core::FE::CellType::line3>(numdofpernode, disname);
    }
    default:
    {
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->shape(), ele->num_node());
      break;
    }
  }

  return nullptr;
}


/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::PoroFluidMultiPhaseEleInterface*
Discret::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::define_problem_type(
    const int numdofpernode, const std::string& disname)
{
  return Discret::ELEMENTS::PoroFluidMultiPhaseEleBoundaryCalc<distype>::instance(
      numdofpernode, disname);
}

FOUR_C_NAMESPACE_CLOSE
