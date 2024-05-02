/*----------------------------------------------------------------------*/
/*! \file

\brief utility class supporting element evaluation for concentrated electrolytes

\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_utils_elch_diffcond.hpp"

#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_newman.hpp"
#include "4C_mat_newman_multiscale.hpp"
#include "4C_scatra_ele_calc_elch_diffcond.hpp"
#include "4C_scatra_ele_calc_elch_diffcond_multiscale.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>*
DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleUtilsElchDiffCond<distype>>(
            new ScaTraEleUtilsElchDiffCond<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::ScaTraEleUtilsElchDiffCond(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelectrode::ScaTraEleUtilsElchElectrode(numdofpernode, numscal, disname)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatElchMat(
    Teuchos::RCP<const CORE::MAT::Material> material, const std::vector<double>& concentrations,
    const double temperature, const INPAR::ELCH::EquPot equpot, const double ffrt,
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanager,
    INPAR::ELCH::DiffCondMat& diffcondmat)
{
  // cast material to electrolyte material
  const Teuchos::RCP<const MAT::ElchMat> elchmat =
      Teuchos::rcp_static_cast<const MAT::ElchMat>(material);

  // safety check
  if (elchmat->NumPhase() != 1)
    FOUR_C_THROW("Can only have a single electrolyte phase at the moment!");

  // extract electrolyte phase
  const Teuchos::RCP<const CORE::MAT::Material> elchphase = elchmat->PhaseById(elchmat->PhaseID(0));

  if (elchphase->MaterialType() == CORE::Materials::m_elchphase)
  {
    // evaluate electrolyte phase
    MatElchPhase(elchphase, concentrations, temperature, equpot, ffrt, diffmanager, diffcondmat);
  }
  else
    FOUR_C_THROW("Invalid material type!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatElchPhase(
    Teuchos::RCP<const CORE::MAT::Material> material, const std::vector<double>& concentrations,
    const double temperature, const INPAR::ELCH::EquPot& equpot, const double& ffrt,
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanager,
    INPAR::ELCH::DiffCondMat& diffcondmat)
{
  // cast material to electrolyte phase
  const Teuchos::RCP<const MAT::ElchPhase> matelchphase =
      Teuchos::rcp_static_cast<const MAT::ElchPhase>(material);

  // set porosity
  diffmanager->SetPhasePoro(matelchphase->Epsilon(), 0);

  // set tortuosity
  diffmanager->SetPhaseTort(matelchphase->Tortuosity(), 0);

  // loop over materials within electrolyte phase
  for (int imat = 0; imat < matelchphase->NumMat(); ++imat)
  {
    const Teuchos::RCP<const CORE::MAT::Material> elchPhaseMaterial =
        matelchphase->MatById(matelchphase->MatID(imat));

    switch (elchPhaseMaterial->MaterialType())
    {
      case CORE::Materials::m_newman:
      case CORE::Materials::m_newman_multiscale:
      {
        // safety check
        if (matelchphase->NumMat() != 1)
          FOUR_C_THROW("Newman material must be the only transported species!");

        // set ion type
        diffcondmat = INPAR::ELCH::diffcondmat_newman;

        // evaluate standard Newman material
        MatNewman(elchPhaseMaterial, concentrations[0], temperature, diffmanager);

        break;
      }

      case CORE::Materials::m_ion:
      {
        // set ion type
        diffcondmat = INPAR::ELCH::diffcondmat_ion;

        myelch::MatIon(elchPhaseMaterial, imat, equpot, diffmanager);

        // calculation of conductivity and transference number based on diffusion coefficient and
        // valence
        if (imat == matelchphase->NumMat() - 1)
        {
          diffmanager->CalcConductivity(matelchphase->NumMat(), ffrt, concentrations);
          diffmanager->CalcTransNum(matelchphase->NumMat(), concentrations);
        }

        break;
      }

      default:
      {
        FOUR_C_THROW("Invalid material type!");
        break;
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::MatNewman(
    Teuchos::RCP<const CORE::MAT::Material> material, const double concentration,
    const double temperature, Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanager)
{
  // cast material to Newman material
  const Teuchos::RCP<const MAT::Newman> matnewman =
      Teuchos::rcp_static_cast<const MAT::Newman>(material);

  // valence of ionic species
  diffmanager->SetValence(matnewman->Valence(), 0);

  // concentration depending diffusion coefficient
  diffmanager->SetIsotropicDiff(
      matnewman->ComputeDiffusionCoefficient(concentration, temperature), 0);
  // derivation of concentration depending diffusion coefficient wrt concentration
  diffmanager->SetConcDerivIsoDiffCoef(
      matnewman->ComputeConcentrationDerivativeOfDiffusionCoefficient(concentration, temperature),
      0, 0);

  // derivation of concentration depending diffusion coefficient wrt temperature
  diffmanager->SetTempDerivIsoDiffCoef(
      matnewman->ComputeTemperatureDerivativeOfDiffusionCoefficient(concentration, temperature), 0,
      0);

  // concentration depending transference number
  diffmanager->SetTransNum(matnewman->ComputeTransferenceNumber(concentration), 0);
  // derivation of concentration depending transference number wrt all ionic species
  diffmanager->SetDerivTransNum(matnewman->ComputeFirstDerivTrans(concentration), 0, 0);

  // thermodynamic factor of electrolyte solution
  diffmanager->SetThermFac(matnewman->ComputeThermFac(concentration));
  // derivative of conductivity with respect to concentrations
  diffmanager->SetDerivThermFac(matnewman->ComputeFirstDerivThermFac(concentration), 0);

  // conductivity and first derivative can maximally depend on one concentration
  // since time curve is used as input routine
  // conductivity of electrolyte solution
  diffmanager->SetCond(matnewman->ComputeConductivity(concentration, temperature));

  // derivative of electronic conductivity w.r.t. concentration
  diffmanager->SetConcDerivCond(
      matnewman->ComputeConcentrationDerivativeOfConductivity(concentration, temperature), 0);

  // derivative of electronic conductivity w.r.t. temperature
  diffmanager->SetTempDerivCond(
      matnewman->ComputeTemperatureDerivativeOfConductivity(concentration, temperature), 0);
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::hex8>;
// template class
// DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::tet10>;
// template class
// DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::pyramid5>;
// template class
// DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
