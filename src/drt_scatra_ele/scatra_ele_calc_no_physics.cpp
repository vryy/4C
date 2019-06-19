/*----------------------------------------------------------------------*/
/*!

\brief Evaluation of a scatra element that does not contain any physics. Currently only implementes
 the minimal set of actions needed for reading the scatra results from a restart file and simulating
 a one-way coupling to the structure. This ImplType is currently not capable to be used in solving
 the scatra equations, as the needed actions are not implemented yet.

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/


#include "scatra_ele_calc_no_physics.h"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalcNoPhysics* delete_me)
{
  static std::map<std::pair<std::string, int>, ScaTraEleCalcNoPhysics<distype, probdim>*> instances;

  std::pair<std::string, int> key(disname, numdofpernode);

  if (delete_me == nullptr)
  {
    if (instances.find(key) == instances.end())
      instances[key] =
          new ScaTraEleCalcNoPhysics<distype, probdim>(numdofpernode, numscal, disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for (typename std::map<std::pair<std::string, int>,
             ScaTraEleCalcNoPhysics<distype, probdim>*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return nullptr;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::Done()
{
  // delete instance
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::ScaTraEleCalcNoPhysics(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs27>;