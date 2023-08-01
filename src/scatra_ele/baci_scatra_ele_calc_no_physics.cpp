/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2


*/
/*----------------------------------------------------------------------*/


#include "baci_scatra_ele_calc_no_physics.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 | singleton access method                                gebauer 06/19 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::pair<std::string, int>>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcNoPhysics<distype, probdim>>(
            new ScaTraEleCalcNoPhysics<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[std::make_pair(disname, numdofpernode)].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
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
#include "baci_scatra_ele_calc_no_physics_fwd.hpp"
