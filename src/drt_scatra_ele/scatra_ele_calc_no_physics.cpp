/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2


*/
/*----------------------------------------------------------------------*/


#include "scatra_ele_calc_no_physics.H"

/*----------------------------------------------------------------------*
 | singleton access method                                gebauer 06/19 |
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
 | singleton destruction                                  gebauer 06/19 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::Done()
{
  // delete instance
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                     gebauer 06/19 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::ScaTraEleCalcNoPhysics(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
}


// include forward declaration of template classes
#include "scatra_ele_calc_no_physics_fwd.hpp"