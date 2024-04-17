/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation for electrodes


\level 2
 */
/*----------------------------------------------------------------------*/
#include "baci_scatra_ele_utils_elch_electrode.hpp"

#include "baci_mat_electrode.hpp"
#include "baci_scatra_ele_calc_elch_electrode.hpp"
#include "baci_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>*
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleUtilsElchElectrode<distype>>(
            new ScaTraEleUtilsElchElectrode<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::ScaTraEleUtilsElchElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelch::ScaTraEleUtilsElch(numdofpernode, numscal, disname)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::MatElectrode(
    Teuchos::RCP<const MAT::Material> material, const double concentration,
    const double temperature, Teuchos::RCP<ScaTraEleDiffManagerElchElectrode> diffmanager)
{
  const auto* matelectrode = static_cast<const MAT::Electrode*>(material.get());

  // diffusion coefficient
  diffmanager->SetIsotropicDiff(
      matelectrode->ComputeDiffusionCoefficient(concentration, temperature), 0);

  // derivative of diffusion coefficient with respect to concentration
  diffmanager->SetConcDerivIsoDiffCoef(
      matelectrode->ComputeConcentrationDerivativeOfDiffusionCoefficient(
          concentration, temperature),
      0, 0);

  // derivative of diffusion coefficient with respect to temperature
  diffmanager->SetTempDerivIsoDiffCoef(
      matelectrode->ComputeTemperatureDerivativeOfDiffusionCoefficient(concentration, temperature),
      0, 0);

  // electronic conductivity
  diffmanager->SetCond(matelectrode->ComputeConductivity(concentration, temperature));

  // derivative of electronic conductivity w.r.t. concentration
  diffmanager->SetConcDerivCond(
      matelectrode->ComputeConcentrationDerivativeOfConductivity(concentration, temperature), 0);

  // derivative of electronic conductivity w.r.t. temperature
  diffmanager->SetTempDerivCond(
      matelectrode->ComputeTemperatureDerivativeOfConductivity(concentration, temperature), 0);
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::hex8>;
// template class
// DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::tet10>;
// template class
// DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::pyramid5>;
// template class
// DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
