/*----------------------------------------------------------------------*/
/*! \file
 \brief factory class providing the implementations of the porofluidmultiphase
        element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/



#include "baci_porofluidmultiphase_ele_factory.H"

#include "baci_global_data.H"
#include "baci_porofluidmultiphase_ele_calc.H"

BACI_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface*
DRT::ELEMENTS::PoroFluidMultiPhaseFactory::ProvideImpl(
    CORE::FE::CellType distype, const int numdofpernode, const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = DRT::Problem::Instance()->NDim();

  switch (distype)
  {
    case CORE::FE::CellType::quad4:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad4>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for quad4 porofluidmultiphase element!");
      break;
    }
      //  case CORE::FE::CellType::quad8:
      //  {
      //    if(ndim==2)
      //      return
      //      DefineProblemType<CORE::FE::CellType::quad8>(numdofpernode,disname);
      //    else
      //      dserror("invalid problem dimension for quad8 porofluidmultiphase element!");
      //    break;
      //  }
    case CORE::FE::CellType::quad9:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad9>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for quad9 porofluidmultiphase element!");
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::tri3>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for tri3 porofluidmultiphase element!");
      break;
    }
      //  case CORE::FE::CellType::tri6:
      //  {
      //    if(ndim==2)
      //      return
      //      DefineProblemType<CORE::FE::CellType::tri6>(numdofpernode,disname);
      //    else
      //      dserror("invalid problem dimension for tri6 porofluidmultiphase element!");
      //    break;
      //  }
    case CORE::FE::CellType::line2:
    {
      if (ndim == 1)
        return DefineProblemType<CORE::FE::CellType::line2>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for line2 porofluidmultiphase element!");
      break;
    }
    case CORE::FE::CellType::line3:
    {
      if (ndim == 1)
        return DefineProblemType<CORE::FE::CellType::line3>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for line3 porofluidmultiphase element!");
      break;
    }
    case CORE::FE::CellType::hex8:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::hex8>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for hex8 porofluidmultiphase element!");
      break;
    }
    case CORE::FE::CellType::hex27:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::hex27>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for hex27 porofluidmultiphase element!");
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tet4>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for tet4 porofluidmultiphase element!");
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tet10>(numdofpernode, disname);
      else
        dserror("invalid problem dimension for tet10 porofluidmultiphase element!");
      break;
    }
    default:
      dserror("Element shape %s not activated. Just do it.",
          CORE::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface*
DRT::ELEMENTS::PoroFluidMultiPhaseFactory::DefineProblemType(
    const int numdofpernode, const std::string& disname)
{
  return DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Instance(numdofpernode, disname);
}

BACI_NAMESPACE_CLOSE
