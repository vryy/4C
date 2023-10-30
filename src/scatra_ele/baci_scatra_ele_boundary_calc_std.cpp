/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#include "baci_scatra_ele_boundary_calc_std.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_discretization_geometry_position_array.H"
#include "baci_lib_globalproblem.H"  // for curves and functions
#include "baci_nurbs_discret.H"
#include "baci_nurbs_discret_nurbs_utils.H"
#include "baci_scatra_ele.H"
#include "baci_scatra_ele_action.H"
#include "baci_scatra_ele_parameter_elch.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcStd<distype, probdim>>(
            new ScaTraEleBoundaryCalcStd<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>::ScaTraEleBoundaryCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<CORE::FE::CellType::nurbs9, 3>;
