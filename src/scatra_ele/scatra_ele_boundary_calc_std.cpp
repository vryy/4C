/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_std.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_action.H"
#include "scatra_ele.H"

#include "lib_globalproblem.H"  // for curves and functions
#include "discretization_fem_general_utils_boundary_integration.H"
#include "discretization_fem_general_utils_fem_shapefunctions.H"
#include "discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "nurbs_discret.H"
#include "nurbs_discret_nurbs_utils.H"
#include "geometry_position_array.H"
#include "utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcStd<distype>>(
            new ScaTraEleBoundaryCalcStd<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      ::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
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
