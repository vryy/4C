/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "4C_elemag_ele_factory.hpp"

#include "4C_elemag_diff_ele_calc.hpp"
#include "4C_elemag_ele_calc.hpp"
#include "4C_elemag_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 |                                                (public) berardocco 02/18 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagEleInterface* Discret::ELEMENTS::ElemagFactory::ProvideImpl(
    Core::FE::CellType distype, std::string problem)
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      return define_problem_type<Core::FE::CellType::hex8>(problem);
    }
    case Core::FE::CellType::hex20:
    {
      return define_problem_type<Core::FE::CellType::hex20>(problem);
    }
    case Core::FE::CellType::hex27:
    {
      return define_problem_type<Core::FE::CellType::hex27>(problem);
    }
    case Core::FE::CellType::tet4:
    {
      return define_problem_type<Core::FE::CellType::tet4>(problem);
    }
    case Core::FE::CellType::tet10:
    {
      return define_problem_type<Core::FE::CellType::tet10>(problem);
    }
    case Core::FE::CellType::wedge6:
    {
      return define_problem_type<Core::FE::CellType::wedge6>(problem);
    }
    /* wedge15 cannot be used since no mesh generator exists
    case Core::FE::CellType::wedge15:
    {
      return define_problem_type<Core::FE::CellType::wedge15>(problem);
    }
    */
    case Core::FE::CellType::pyramid5:
    {
      return define_problem_type<Core::FE::CellType::pyramid5>(problem);
    }
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
    // Nurbs support
    case Core::FE::CellType::nurbs9:
    {
      return define_problem_type<Core::FE::CellType::nurbs9>(problem);
    }
    case Core::FE::CellType::nurbs27:
    {
      return define_problem_type<Core::FE::CellType::nurbs27>(problem);
    }
    // no 1D elements
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          Core::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                (public) berardocco 02/18 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ElemagEleInterface* Discret::ELEMENTS::ElemagFactory::define_problem_type(
    std::string problem)
{
  if (problem == "std")
    return Discret::ELEMENTS::ElemagEleCalc<distype>::Instance();
  else if (problem == "diff")
    return Discret::ELEMENTS::ElemagDiffEleCalc<distype>::Instance();
  else
    FOUR_C_THROW("Defined problem type does not exist!!");

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
