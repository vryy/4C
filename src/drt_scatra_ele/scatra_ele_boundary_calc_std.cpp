/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_std.cpp

\brief evaluation of scatra boundary terms at integration points

\level 1

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
</pre>
 */
/*----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_std.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_action.H"
#include "scatra_ele.H"

#include "../drt_lib/drt_globalproblem.H"  // for curves and functions
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleBoundaryCalcStd* delete_me)
{
  static std::map<std::string, ScaTraEleBoundaryCalcStd<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcStd<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleBoundaryCalcStd<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::ScaTraEleBoundaryCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<DRT::Element::nurbs9>;
