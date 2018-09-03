/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_std.cpp

\brief evaluation of scalar transport elements for standard scalar transport problems

\level 1

<pre>
\maintainer Rui Fang
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
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalcStd* delete_me)
{
  static std::map<std::pair<std::string, int>, ScaTraEleCalcStd<distype, probdim>*> instances;

  std::pair<std::string, int> key(disname, numdofpernode);

  if (delete_me == NULL)
  {
    if (instances.find(key) == instances.end())
      instances[key] = new ScaTraEleCalcStd<distype, probdim>(numdofpernode, numscal, disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for (typename std::map<std::pair<std::string, int>,
             ScaTraEleCalcStd<distype, probdim>*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::Done()
{
  // delete instance
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::ScaTraEleCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs27>;
