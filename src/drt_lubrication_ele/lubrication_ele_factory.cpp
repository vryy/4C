/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_factory.cpp

\brief

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "../drt_lubrication_ele/lubrication_ele_factory.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lubrication_ele/lubrication_ele_calc.H"

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
    case DRT::Element::quad4:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::quad4, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::quad4, 3>(disname);
      else
        dserror("invalid problem dimension for quad4 lubrication element!");
      break;
    }
    case DRT::Element::quad8:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::quad8, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::quad8, 3>(disname);
      else
        dserror("invalid problem dimension for quad8 lubrication element!");
      break;
    }
    case DRT::Element::quad9:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::quad9, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::quad9, 3>(disname);
      else
        dserror("invalid problem dimension for quad9 lubrication element!");
      break;
    }
    case DRT::Element::tri3:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::tri3, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::tri3, 3>(disname);
      else
        dserror("invalid problem dimension for tri3 lubrication element!");
      break;
    }
    case DRT::Element::tri6:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::tri6, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::tri6, 3>(disname);
      else
        dserror("invalid problem dimension for tri6 lubrication element!");
      break;
    }
    case DRT::Element::line2:
    {
      if (ndim == 1)
        return DefineProblemType<DRT::Element::line2, 1>(disname);
      else if (ndim == 2)
        return DefineProblemType<DRT::Element::line2, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::line2, 3>(disname);
      else
        dserror("invalid problem dimension for line2 lubrication element!");
      break;
    }
    case DRT::Element::line3:
    {
      if (ndim == 1)
        return DefineProblemType<DRT::Element::line3, 1>(disname);
      else if (ndim == 2)
        return DefineProblemType<DRT::Element::line3, 2>(disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::line3, 3>(disname);
      else
        dserror("invalid problem dimension for line3 lubrication element!");
      break;
    }
    default:
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
  }
  return NULL;
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
