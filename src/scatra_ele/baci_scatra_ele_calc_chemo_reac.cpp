/*----------------------------------------------------------------------*/
/*! \file
 \brief main file containing routines for calculation of scatra element with chemotactic AND
reactive scalars

\level 2

 *----------------------------------------------------------------------*/

#include "baci_scatra_ele_calc_chemo_reac.hpp"

#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_element.hpp"
#include "baci_mat_list.hpp"
#include "baci_mat_list_chemoreac.hpp"
#include "baci_mat_list_chemotaxis.hpp"
#include "baci_mat_list_reactions.hpp"
#include "baci_mat_scatra.hpp"
#include "baci_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcChemoReac<distype, probdim>::ScaTraEleCalcChemoReac(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcChemo<distype, probdim>::ScaTraEleCalcChemo(
          numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype, probdim>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcChemoReac<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcChemoReac<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcChemoReac<distype, probdim>>(
            new ScaTraEleCalcChemoReac<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                    thon 06/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcChemoReac<distype, probdim>::GetMaterialParams(
    const DRT::Element* ele,      //!< the element we are dealing with
    std::vector<double>& densn,   //!< density at t_(n)
    std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,  //!< density at t_(n+alpha_M)
    double& visc,                 //!< fluid viscosity
    const int iquad               //!< id of current gauss point (default = -1)
)
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // We may have some reactive and some non-reactive elements in one discretisation.
  // But since the calculation classes are singleton, we have to reset all reactive stuff in case
  // of non-reactive elements:
  advreac::ReaManager()->Clear(my::numscal_);

  // We may have some chemotactic and some non-chemotactic discretisation.
  // But since the calculation classes are singleton, we have to reset all chemotaxis stuff each
  // time
  chemo::ClearChemotaxisTerms();

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList> actmat =
        Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0; k < my::numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      my::Materials(singlemat, k, densn[k], densnp[k], densam[k], visc, iquad);
    }
  }

  else if (material->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    const Teuchos::RCP<MAT::MatListReactions> actmat =
        Teuchos::rcp_dynamic_cast<MAT::MatListReactions>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0; k < my::numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      // Note: order is important here!!
      advreac::Materials(singlemat, k, densn[k], densnp[k], densam[k], visc, iquad);

      advreac::SetAdvancedReactionTerms(
          k, actmat, advreac::GetGpCoord());  // every reaction calculation stuff happens in here!!
    }
  }

  else if (material->MaterialType() == INPAR::MAT::m_matlist_chemotaxis)
  {
    const Teuchos::RCP<MAT::MatListChemotaxis> actmat =
        Teuchos::rcp_dynamic_cast<MAT::MatListChemotaxis>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    chemo::GetChemotaxisCoefficients(
        material);  // read all chemotaxis input from material and copy it into local variables

    for (int k = 0; k < my::numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      my::Materials(singlemat, k, densn[k], densnp[k], densam[k], visc, iquad);
    }
  }

  else if (material->MaterialType() == INPAR::MAT::m_matlist_chemoreac)
  {
    const Teuchos::RCP<MAT::MatListReactions> actmat =
        Teuchos::rcp_dynamic_cast<MAT::MatListReactions>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    chemo::GetChemotaxisCoefficients(
        material);  // read all chemotaxis input from material and copy it into local variables

    for (int k = 0; k < my::numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      // Note: order is important here!!
      my::Materials(singlemat, k, densn[k], densnp[k], densam[k], visc, iquad);
      advreac::SetAdvancedReactionTerms(
          k, actmat, advreac::GetGpCoord());  // every reaction calculation stuff happens in here!!
    }
  }

  else
  {
    advreac::Materials(material, 0, densn[0], densnp[0], densam[0], visc, iquad);
  }

  return;
}  // ScaTraEleCalc::GetMaterialParams


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::nurbs9>;
// template class DRT::ELEMENTS::ScaTraEleCalcChemoReac<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
