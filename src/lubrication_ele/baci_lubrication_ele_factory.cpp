/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "baci_lubrication_ele_factory.H"

#include "baci_lib_globalproblem.H"
#include "baci_lubrication_ele_calc.H"

/*--------------------------------------------------------------------------*
 |                                                     (public) wirtz 10/15 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::LubricationEleInterface* DRT::ELEMENTS::LubricationFactory::ProvideImpl(
    DRT::Element::DiscretizationType distype, const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = DRT::Problem::Instance()->NDim();

  switch (distype)
  {
    case DRT::Element::DiscretizationType::quad4:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::quad4, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::quad4, 3>(disname);
      else
        dserror("invalid problem dimension for quad4 lubrication element!");
      break;
    }
    case DRT::Element::DiscretizationType::quad8:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::quad8, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::quad8, 3>(disname);
      else
        dserror("invalid problem dimension for quad8 lubrication element!");
      break;
    }
    case DRT::Element::DiscretizationType::quad9:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::quad9, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::quad9, 3>(disname);
      else
        dserror("invalid problem dimension for quad9 lubrication element!");
      break;
    }
    case DRT::Element::DiscretizationType::tri3:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::tri3, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::tri3, 3>(disname);
      else
        dserror("invalid problem dimension for tri3 lubrication element!");
      break;
    }
    case DRT::Element::DiscretizationType::tri6:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::tri6, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::tri6, 3>(disname);
      else
        dserror("invalid problem dimension for tri6 lubrication element!");
      break;
    }
    case DRT::Element::DiscretizationType::line2:
    {
      if (ndim == 1)
        return DefineProblemType<DRT::Element::DiscretizationType::line2, 1>(disname);
      else if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::line2, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::line2, 3>(disname);
      else
        dserror("invalid problem dimension for line2 lubrication element!");
      break;
    }
    case DRT::Element::DiscretizationType::line3:
    {
      if (ndim == 1)
        return DefineProblemType<DRT::Element::DiscretizationType::line3, 1>(disname);
      else if (ndim == 2)
        return DefineProblemType<DRT::Element::DiscretizationType::line3, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::DiscretizationType::line3, 3>(disname);
      else
        dserror("invalid problem dimension for line3 lubrication element!");
      break;
    }
    default:
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                     (public) wirtz 10/15 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::LubricationEleInterface* DRT::ELEMENTS::LubricationFactory::DefineProblemType(
    const std::string& disname)
{
  return DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Instance(disname);
}
