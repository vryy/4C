/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_std.cpp

\brief evaluation of scalar transport elements for standard scalar transport problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_std.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcStd<distype>* DRT::ELEMENTS::ScaTraEleCalcStd<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    bool create
    )
{
  static ScaTraEleCalcStd<distype>* instance;

  if(create)
  {
    if(instance == NULL)
      instance = new ScaTraEleCalcStd<distype>(numdofpernode,numscal);
  }

  else if(instance != NULL)
  {
    delete instance;
    instance = NULL;
  }

  return instance;
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcStd<distype>::Done()
{
  // delete instance
  Instance(0,0,false);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcStd<distype>::ScaTraEleCalcStd(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{
  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs27>;
