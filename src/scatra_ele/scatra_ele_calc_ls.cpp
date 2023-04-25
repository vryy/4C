/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluations for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_ls.H"
#include "scatra_ele_parameter_std.H"
#include "headers_singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>* DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcLS<distype>>(
            new ScaTraEleCalcLS<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      ::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>::ScaTraEleCalcLS(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // safety check
  if (my::scatrapara_->RBSubGrVel()) dserror("CalcSubgrVelocityLevelSet not available anymore");

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs27>;
