/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_calc_std.cpp

\brief standard routines for calculation of fluid element

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*/
/*----------------------------------------------------------------------*/

#include "meshfree_fluid_cell_calc_std.H"
#include "../drt_fluid_ele/fluid_ele_parameter_std.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidCellCalcStd<distype> * DRT::ELEMENTS::MeshfreeFluidCellCalcStd<distype>::Instance( bool create )
{
  static MeshfreeFluidCellCalcStd<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new MeshfreeFluidCellCalcStd<distype>();
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
void DRT::ELEMENTS::MeshfreeFluidCellCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidCellCalcStd<distype>::MeshfreeFluidCellCalcStd()
  : DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::MeshfreeFluidCellCalc()
{
  my::fldpara_=DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

// template classes
template class DRT::ELEMENTS::MeshfreeFluidCellCalcStd<DRT::Element::hex8>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalcStd<DRT::Element::tet4>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalcStd<DRT::Element::quad4>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalcStd<DRT::Element::tri3>;

