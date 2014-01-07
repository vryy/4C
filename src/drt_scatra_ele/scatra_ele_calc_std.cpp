/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_std.cpp

\brief Factory of scatra elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_std.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcStd<distype> * DRT::ELEMENTS::ScaTraEleCalcStd<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcStd<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcStd<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcStd<distype>::ScaTraEleCalcStd(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{

}


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad8>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs27>;



