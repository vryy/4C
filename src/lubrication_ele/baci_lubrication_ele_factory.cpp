/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "baci_lubrication_ele_factory.H"

#include "baci_lib_globalproblem.H"
#include "baci_lubrication_ele_calc.H"

BACI_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 |                                                     (public) wirtz 10/15 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::LubricationEleInterface* DRT::ELEMENTS::LubricationFactory::ProvideImpl(
    CORE::FE::CellType distype, const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = DRT::Problem::Instance()->NDim();

  switch (distype)
  {
    case CORE::FE::CellType::quad4:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad4, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::quad4, 3>(disname);
      else
        dserror("invalid problem dimension for quad4 lubrication element!");
      break;
    }
    case CORE::FE::CellType::quad8:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad8, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::quad8, 3>(disname);
      else
        dserror("invalid problem dimension for quad8 lubrication element!");
      break;
    }
    case CORE::FE::CellType::quad9:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad9, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::quad9, 3>(disname);
      else
        dserror("invalid problem dimension for quad9 lubrication element!");
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::tri3, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tri3, 3>(disname);
      else
        dserror("invalid problem dimension for tri3 lubrication element!");
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::tri6, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tri6, 3>(disname);
      else
        dserror("invalid problem dimension for tri6 lubrication element!");
      break;
    }
    case CORE::FE::CellType::line2:
    {
      if (ndim == 1)
        return DefineProblemType<CORE::FE::CellType::line2, 1>(disname);
      else if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::line2, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::line2, 3>(disname);
      else
        dserror("invalid problem dimension for line2 lubrication element!");
      break;
    }
    case CORE::FE::CellType::line3:
    {
      if (ndim == 1)
        return DefineProblemType<CORE::FE::CellType::line3, 1>(disname);
      else if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::line3, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::line3, 3>(disname);
      else
        dserror("invalid problem dimension for line3 lubrication element!");
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
 |                                                     (public) wirtz 10/15 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::LubricationEleInterface* DRT::ELEMENTS::LubricationFactory::DefineProblemType(
    const std::string& disname)
{
  return DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Instance(disname);
}

BACI_NAMESPACE_CLOSE
