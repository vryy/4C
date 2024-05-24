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
DRT::ELEMENTS::ElemagEleInterface* DRT::ELEMENTS::ElemagFactory::ProvideImpl(
    CORE::FE::CellType distype, std::string problem)
{
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      return define_problem_type<CORE::FE::CellType::hex8>(problem);
    }
    case CORE::FE::CellType::hex20:
    {
      return define_problem_type<CORE::FE::CellType::hex20>(problem);
    }
    case CORE::FE::CellType::hex27:
    {
      return define_problem_type<CORE::FE::CellType::hex27>(problem);
    }
    case CORE::FE::CellType::tet4:
    {
      return define_problem_type<CORE::FE::CellType::tet4>(problem);
    }
    case CORE::FE::CellType::tet10:
    {
      return define_problem_type<CORE::FE::CellType::tet10>(problem);
    }
    case CORE::FE::CellType::wedge6:
    {
      return define_problem_type<CORE::FE::CellType::wedge6>(problem);
    }
    /* wedge15 cannot be used since no mesh generator exists
    case CORE::FE::CellType::wedge15:
    {
      return define_problem_type<CORE::FE::CellType::wedge15>(problem);
    }
    */
    case CORE::FE::CellType::pyramid5:
    {
      return define_problem_type<CORE::FE::CellType::pyramid5>(problem);
    }
    case CORE::FE::CellType::quad4:
    {
      return define_problem_type<CORE::FE::CellType::quad4>(problem);
    }
    case CORE::FE::CellType::quad8:
    {
      return define_problem_type<CORE::FE::CellType::quad8>(problem);
    }
    case CORE::FE::CellType::quad9:
    {
      return define_problem_type<CORE::FE::CellType::quad9>(problem);
    }
    case CORE::FE::CellType::tri3:
    {
      return define_problem_type<CORE::FE::CellType::tri3>(problem);
    }
    case CORE::FE::CellType::tri6:
    {
      return define_problem_type<CORE::FE::CellType::tri6>(problem);
    }
    // Nurbs support
    case CORE::FE::CellType::nurbs9:
    {
      return define_problem_type<CORE::FE::CellType::nurbs9>(problem);
    }
    case CORE::FE::CellType::nurbs27:
    {
      return define_problem_type<CORE::FE::CellType::nurbs27>(problem);
    }
    // no 1D elements
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          CORE::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                (public) berardocco 02/18 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagEleInterface* DRT::ELEMENTS::ElemagFactory::define_problem_type(
    std::string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::ElemagEleCalc<distype>::Instance();
  else if (problem == "diff")
    return DRT::ELEMENTS::ElemagDiffEleCalc<distype>::Instance();
  else
    FOUR_C_THROW("Defined problem type does not exist!!");

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
