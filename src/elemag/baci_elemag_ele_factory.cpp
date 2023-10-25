/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "baci_elemag_ele_factory.H"

#include "baci_elemag_diff_ele_calc.H"
#include "baci_elemag_ele_calc.H"
#include "baci_elemag_ele_interface.H"

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
      return DefineProblemType<CORE::FE::CellType::hex8>(problem);
    }
    case CORE::FE::CellType::hex20:
    {
      return DefineProblemType<CORE::FE::CellType::hex20>(problem);
    }
    case CORE::FE::CellType::hex27:
    {
      return DefineProblemType<CORE::FE::CellType::hex27>(problem);
    }
    case CORE::FE::CellType::tet4:
    {
      return DefineProblemType<CORE::FE::CellType::tet4>(problem);
    }
    case CORE::FE::CellType::tet10:
    {
      return DefineProblemType<CORE::FE::CellType::tet10>(problem);
    }
    case CORE::FE::CellType::wedge6:
    {
      return DefineProblemType<CORE::FE::CellType::wedge6>(problem);
    }
    /* wedge15 cannot be used since no mesh generator exists
    case CORE::FE::CellType::wedge15:
    {
      return DefineProblemType<CORE::FE::CellType::wedge15>(problem);
    }
    */
    case CORE::FE::CellType::pyramid5:
    {
      return DefineProblemType<CORE::FE::CellType::pyramid5>(problem);
    }
    case CORE::FE::CellType::quad4:
    {
      return DefineProblemType<CORE::FE::CellType::quad4>(problem);
    }
    case CORE::FE::CellType::quad8:
    {
      return DefineProblemType<CORE::FE::CellType::quad8>(problem);
    }
    case CORE::FE::CellType::quad9:
    {
      return DefineProblemType<CORE::FE::CellType::quad9>(problem);
    }
    case CORE::FE::CellType::tri3:
    {
      return DefineProblemType<CORE::FE::CellType::tri3>(problem);
    }
    case CORE::FE::CellType::tri6:
    {
      return DefineProblemType<CORE::FE::CellType::tri6>(problem);
    }
    // Nurbs support
    case CORE::FE::CellType::nurbs9:
    {
      return DefineProblemType<CORE::FE::CellType::nurbs9>(problem);
    }
    case CORE::FE::CellType::nurbs27:
    {
      return DefineProblemType<CORE::FE::CellType::nurbs27>(problem);
    }
    // no 1D elements
    default:
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                (public) berardocco 02/18 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagEleInterface* DRT::ELEMENTS::ElemagFactory::DefineProblemType(
    std::string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::ElemagEleCalc<distype>::Instance();
  else if (problem == "diff")
    return DRT::ELEMENTS::ElemagDiffEleCalc<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return nullptr;
}
