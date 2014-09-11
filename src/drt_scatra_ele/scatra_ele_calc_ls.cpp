/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_ls.cpp

\brief evaluations for level sets

<pre>
Maintainer: Ursula Rasthofer
            erasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_ls.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype> * DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcLS<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcLS<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>::ScaTraEleCalcLS(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{

}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs27>;



