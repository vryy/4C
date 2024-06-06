/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation for electrodes


\level 2
 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_utils_elch_electrode.hpp"

#include "4C_mat_electrode.hpp"
#include "4C_scatra_ele_calc_elch_electrode.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>*
Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleUtilsElchElectrode<distype>>(
            new ScaTraEleUtilsElchElectrode<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::ScaTraEleUtilsElchElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelch::ScaTraEleUtilsElch(numdofpernode, numscal, disname)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::mat_electrode(
    Teuchos::RCP<const Core::Mat::Material> material, const double concentration,
    const double temperature, Teuchos::RCP<ScaTraEleDiffManagerElchElectrode> diffmanager)
{
  const auto* matelectrode = static_cast<const Mat::Electrode*>(material.get());

  // diffusion coefficient
  diffmanager->SetIsotropicDiff(
      matelectrode->compute_diffusion_coefficient(concentration, temperature), 0);

  // derivative of diffusion coefficient with respect to concentration
  diffmanager->set_conc_deriv_iso_diff_coef(
      matelectrode->compute_concentration_derivative_of_diffusion_coefficient(
          concentration, temperature),
      0, 0);

  // derivative of diffusion coefficient with respect to temperature
  diffmanager->set_temp_deriv_iso_diff_coef(
      matelectrode->compute_temperature_derivative_of_diffusion_coefficient(
          concentration, temperature),
      0, 0);

  // electronic conductivity
  diffmanager->SetCond(matelectrode->compute_conductivity(concentration, temperature));

  // derivative of electronic conductivity w.r.t. concentration
  diffmanager->SetConcDerivCond(
      matelectrode->compute_concentration_derivative_of_conductivity(concentration, temperature),
      0);

  // derivative of electronic conductivity w.r.t. temperature
  diffmanager->SetTempDerivCond(
      matelectrode->compute_temperature_derivative_of_conductivity(concentration, temperature), 0);
}


// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::nurbs3>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::hex8>;
// template class
// Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::tet10>;
// template class
// Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::pyramid5>;
// template class
// Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
