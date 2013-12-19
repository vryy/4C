/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.cpp

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_elch.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElch<distype> * DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcElch<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcElch<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElch<distype>::ScaTraEleCalcElch(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{

}



// template classes
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
////template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
////template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
////template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad9>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri3>;
////template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;



