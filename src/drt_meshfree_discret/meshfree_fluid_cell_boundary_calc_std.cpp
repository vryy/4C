/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_boundary_calc_std.cpp

\brief standard routines for calculation of meshfree fluid boundary cell

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#include "meshfree_fluid_cell_boundary_calc_std.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_fluid_ele/fluid_ele_parameter_std.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<distype>* DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<distype>::Instance( bool create )
{
  static MeshfreeFluidBoundaryCalcStd<distype>* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new MeshfreeFluidBoundaryCalcStd<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<distype>::MeshfreeFluidBoundaryCalcStd()
  : DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<distype>::MeshfreeFluidBoundaryCalc()
{
  // pointer to class FluidImplParameter
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<DRT::Element::quad4>;
template class DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<DRT::Element::tri3>;
template class DRT::ELEMENTS::MeshfreeFluidBoundaryCalcStd<DRT::Element::line2>;
