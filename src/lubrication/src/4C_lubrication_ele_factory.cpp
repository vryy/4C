// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_lubrication_ele_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_lubrication_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 |                                                     (public) wirtz 10/15 |
 *--------------------------------------------------------------------------*/
Discret::Elements::LubricationEleInterface* Discret::Elements::LubricationFactory::provide_impl(
    Core::FE::CellType distype, const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = Global::Problem::instance()->n_dim();

  switch (distype)
  {
    case Core::FE::CellType::quad4:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::quad4, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::quad4, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for quad4 lubrication element!");
      break;
    }
    case Core::FE::CellType::quad8:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::quad8, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::quad8, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for quad8 lubrication element!");
      break;
    }
    case Core::FE::CellType::quad9:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::quad9, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::quad9, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for quad9 lubrication element!");
      break;
    }
    case Core::FE::CellType::tri3:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::tri3, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::tri3, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for tri3 lubrication element!");
      break;
    }
    case Core::FE::CellType::tri6:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::tri6, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::tri6, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for tri6 lubrication element!");
      break;
    }
    case Core::FE::CellType::line2:
    {
      if (ndim == 1)
        return define_problem_type<Core::FE::CellType::line2, 1>(disname);
      else if (ndim == 2)
        return define_problem_type<Core::FE::CellType::line2, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::line2, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for line2 lubrication element!");
      break;
    }
    case Core::FE::CellType::line3:
    {
      if (ndim == 1)
        return define_problem_type<Core::FE::CellType::line3, 1>(disname);
      else if (ndim == 2)
        return define_problem_type<Core::FE::CellType::line3, 2>(disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::line3, 3>(disname);
      else
        FOUR_C_THROW("invalid problem dimension for line3 lubrication element!");
      break;
    }
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          Core::FE::cell_type_to_string(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                     (public) wirtz 10/15 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::Elements::LubricationEleInterface*
Discret::Elements::LubricationFactory::define_problem_type(const std::string& disname)
{
  return Discret::Elements::LubricationEleCalc<distype, probdim>::instance(disname);
}

FOUR_C_NAMESPACE_CLOSE
