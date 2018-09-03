/*!----------------------------------------------------------------------
\file fluid_ele_boundary_calc.H

\brief evaluation of fluid terms at integration points

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*----------------------------------------------------------------------*/

#include "fluid_ele_boundary_calc_std.H"
//#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementtype.H"
//#include "fluid_ele.H"
//#include "fluid_ele_action.H"
#include "fluid_ele_parameter_std.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>*
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::Instance(bool create)
{
  static FluidEleBoundaryCalcStd<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new FluidEleBoundaryCalcStd<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::FluidEleBoundaryCalcStd()
    : DRT::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
{
  // pointer to class FluidImplParameter
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::line2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::line3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::nurbs2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::nurbs4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::nurbs9>;
