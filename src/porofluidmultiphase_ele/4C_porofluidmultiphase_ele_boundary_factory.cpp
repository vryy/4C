/*----------------------------------------------------------------------*/
/*! \file
 \brief factory class providing the implementations of the porofluidmultiphase
        boundary element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_porofluidmultiphase_ele_boundary_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_element.hpp"
#include "4C_porofluidmultiphase_ele_boundary_calc.hpp"
#include "4C_porofluidmultiphase_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN


/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface*
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::ProvideImpl(
    const DRT::Element* ele, const int numdofpernode, const std::string& disname)
{
  switch (ele->Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return DefineProblemType<CORE::FE::CellType::quad4>(numdofpernode, disname);
    }
    case CORE::FE::CellType::quad8:
    {
      return DefineProblemType<CORE::FE::CellType::quad8>(numdofpernode, disname);
    }
    case CORE::FE::CellType::quad9:
    {
      return DefineProblemType<CORE::FE::CellType::quad9>(numdofpernode, disname);
    }
    case CORE::FE::CellType::tri3:
    {
      return DefineProblemType<CORE::FE::CellType::tri3>(numdofpernode, disname);
    }
    case CORE::FE::CellType::tri6:
    {
      return DefineProblemType<CORE::FE::CellType::tri6>(numdofpernode, disname);
    }
    case CORE::FE::CellType::line2:
    {
      return DefineProblemType<CORE::FE::CellType::line2>(numdofpernode, disname);
    }
    case CORE::FE::CellType::line3:
    {
      return DefineProblemType<CORE::FE::CellType::line3>(numdofpernode, disname);
    }
    default:
    {
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
      break;
    }
  }

  return nullptr;
}


/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface*
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::DefineProblemType(
    const int numdofpernode, const std::string& disname)
{
  return DRT::ELEMENTS::PoroFluidMultiPhaseEleBoundaryCalc<distype>::Instance(
      numdofpernode, disname);
}

FOUR_C_NAMESPACE_CLOSE
