/*!----------------------------------------------------------------------
\file fluid_ele_calc_immersed.cpp

\brief calc class for immersed problems

<pre>
Maintainers: Andreas Rauch & Anh-Tu Vuong
             {rauch,vuong}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240 / 15264
</pre>
*----------------------------------------------------------------------*/

#include "fluid_ele.H"
#include "fluid_ele_utils.H"
#include "fluid_ele_action.H"
#include "fluid_ele_calc.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_calc_immersed.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcImmersed<distype> * DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Instance( bool create )
{
  static FluidEleCalcImmersed<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcImmersed<distype>();
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
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcImmersed<distype>::FluidEleCalcImmersed()
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
{
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Evaluate(
    DRT::ELEMENTS::Fluid*                       ele,
    DRT::Discretization &                       discretization,
    const std::vector<int> &                    lm,
    Teuchos::ParameterList&                     params,
    Teuchos::RCP<MAT::Material> &               mat,
    Epetra_SerialDenseMatrix&                   elemat1_epetra,
    Epetra_SerialDenseMatrix&                   elemat2_epetra,
    Epetra_SerialDenseVector&                   elevec1_epetra,
    Epetra_SerialDenseVector&                   elevec2_epetra,
    Epetra_SerialDenseVector&                   elevec3_epetra,
    const DRT::UTILS::GaussIntegration &        intpoints)
{
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::nurbs27>;
