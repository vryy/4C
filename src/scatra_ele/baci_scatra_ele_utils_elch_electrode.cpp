/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation for electrodes


\level 2
 */
/*----------------------------------------------------------------------*/
#include "baci_scatra_ele_utils_elch_electrode.H"
#include "baci_scatra_ele_calc_elch_electrode.H"

#include "baci_mat_electrode.H"
#include "baci_utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::ScaTraEleUtilsElchElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelch::ScaTraEleUtilsElch(numdofpernode, numscal, disname)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
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
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<DRT::Element::nurbs27>;
