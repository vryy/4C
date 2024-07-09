/*----------------------------------------------------------------------*/
/*! \file

\brief factory class into templated evaluators for fluid boundary integration

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_boundary_factory.hpp"

#include "4C_fluid_ele_boundary_calc_poro.hpp"
#include "4C_fluid_ele_boundary_calc_std.hpp"
#include "4C_fluid_ele_boundary_interface.hpp"
#include "4C_fluid_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::FluidBoundaryInterface* Discret::ELEMENTS::FluidBoundaryFactory::provide_impl(
    Core::FE::CellType distype, std::string problem)
{
  switch (distype)
  {
    case Core::FE::CellType::quad4:
    {
      return define_problem_type<Core::FE::CellType::quad4>(problem);
    }
    case Core::FE::CellType::quad8:
    {
      return define_problem_type<Core::FE::CellType::quad8>(problem);
    }
    case Core::FE::CellType::quad9:
    {
      return define_problem_type<Core::FE::CellType::quad9>(problem);
    }
    case Core::FE::CellType::tri3:
    {
      return define_problem_type<Core::FE::CellType::tri3>(problem);
    }
    case Core::FE::CellType::tri6:
    {
      return define_problem_type<Core::FE::CellType::tri6>(problem);
    }
    case Core::FE::CellType::line2:
    {
      return define_problem_type<Core::FE::CellType::line2>(problem);
    }
    case Core::FE::CellType::line3:
    {
      return define_problem_type<Core::FE::CellType::line3>(problem);
    }
    case Core::FE::CellType::nurbs2:
    {
      return define_problem_type<Core::FE::CellType::nurbs2>(problem);
    }
    case Core::FE::CellType::nurbs3:
    {
      return define_problem_type<Core::FE::CellType::nurbs3>(problem);
    }
    case Core::FE::CellType::nurbs4:
    {
      return define_problem_type<Core::FE::CellType::nurbs4>(problem);
    }
    case Core::FE::CellType::nurbs9:
    {
      return define_problem_type<Core::FE::CellType::nurbs9>(problem);
    }
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          Core::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidBoundaryInterface*
Discret::ELEMENTS::FluidBoundaryFactory::define_problem_type(std::string problem)
{
  if (problem == "std")
    return Discret::ELEMENTS::FluidEleBoundaryCalcStd<distype>::instance();
  else if (problem == "poro")
    return Discret::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::instance();
  else if (problem == "poro_p1")
    return Discret::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>::instance();
  else
    FOUR_C_THROW("Defined problem type does not exist!!");

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
