/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#include "4C_scatra_ele_boundary_calc_std.hpp"

#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_global_data.hpp"  // for curves and functions
#include "4C_nurbs_discret.hpp"
#include "4C_nurbs_discret_nurbs_utils.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>*
Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcStd<distype, probdim>>(
            new ScaTraEleBoundaryCalcStd<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>::ScaTraEleBoundaryCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::quad4, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::quad8, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::quad9, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::tri6, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::line3, 2>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::nurbs3, 2>;
template class Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<Core::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
