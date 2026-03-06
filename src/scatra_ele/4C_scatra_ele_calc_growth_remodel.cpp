// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_calc_growth_remodel.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_mixture.hpp"
#include "4C_mat_scatra_growth_remodel.hpp"
#include "4C_mixture_constituent_remodelfiber_ssi.hpp"
#include "4C_scatra_ele_calc.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::ScaTraEleCalcGrowthRemodel(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::Elements::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname)
{
}

template <Core::FE::CellType distype, int probdim>
Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>*
Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::pair<std::string, int>>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcGrowthRemodel<distype, probdim>>(
            new ScaTraEleCalcGrowthRemodel<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[std::make_pair(disname, numdofpernode)].instance(
      Core::Utils::SingletonAction::create, numdofpernode, numscal, disname);
}

template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::get_material_params(
    const Core::Elements::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // get the material
  std::shared_ptr<Core::Mat::Material> material = ele->material();

  rea_manager()->clear(my::numscal_);

  if (material->material_type() == Core::Materials::m_matlist)
  {
    const std::shared_ptr<Mat::MatList> actmat = std::dynamic_pointer_cast<Mat::MatList>(material);
    if (actmat->num_mat() < my::numscal_) FOUR_C_THROW("Not enough materials in MatList.");

    for (int k = 0; k < my::numscal_; ++k)
    {
      int matid = actmat->mat_id(k);
      const std::shared_ptr<Core::Mat::Material> singlemat = actmat->material_by_id(matid);
      materials(singlemat, k, densn[k], densnp[k], densam[k], visc, iquad);
    }
  }
  else
  {
    materials(material, 0, densn[0], densnp[0], densam[0], visc, iquad);
  }
}

template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::materials(
    const std::shared_ptr<const Core::Mat::Material> material, const int k, double& densn,
    double& densnp, double& densam, double& visc, const int iquad)
{
  switch (material->material_type())
  {
    case Core::Materials::m_scatra_gr:
    {
      mat_scatra_gr(material, k, densn, densnp, densam, visc, iquad);
      break;
    }
    default:
    {
      my::materials(material, k, densn, densnp, densam, visc, iquad);
      break;
    }
  }
}

template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::mat_scatra_gr(
    const std::shared_ptr<const Core::Mat::Material> material, const int k, double& densn,
    double& densnp, double& densam, double& visc, const int iquad)
{
  const auto actmat = std::dynamic_pointer_cast<const Mat::ScatraGrowthRemodelMat>(material);

  my::diffmanager_->set_isotropic_diff(actmat->diffusivity(), k);

  FOUR_C_ASSERT(
      constituent_[k] != nullptr, "constituent_[{}] is null. Was setup_calc() called?", k);

  const double rea_coeff =
      get_stress_dependent_reaction_coeff(*constituent_[k], actmat->scalar_quantity(), iquad);
  my::reamanager_->set_rea_coeff(rea_coeff, k);
}

template <Core::FE::CellType distype, int probdim>
int Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::setup_calc(
    Core::Elements::Element* ele, Core::FE::Discretization& discretization)
{
  // call base class setup
  const int result = my::setup_calc(ele, discretization);
  if (result == -1) return result;

  // resolve structural material once per element
  const std::shared_ptr<Mat::Mixture> structmat = get_struct_material(ele);

  constituent_.assign(my::numscal_, nullptr);

  const std::shared_ptr<Core::Mat::Material> scatra_material = ele->material(0);

  auto resolve = [&](const std::shared_ptr<const Core::Mat::Material>& mat)
      -> const Mixture::MixtureConstituentRemodelFiberSsi*
  {
    if (mat->material_type() != Core::Materials::m_scatra_gr) return nullptr;

    const auto actmat = std::dynamic_pointer_cast<const Mat::ScatraGrowthRemodelMat>(mat);
    const int struct_mat_id = actmat->structure_material_id();

    const auto* constituent = structmat->ssi_constituent_by_material_id(struct_mat_id);

    FOUR_C_ASSERT_ALWAYS(constituent != nullptr,
        "Constituent with material id {} is not a MixtureConstituentRemodelFiberSsi. "
        "ScaTraEleCalcGrowthRemodel requires SSI RemodelFiber constituents for m_scatra_gr "
        "materials.",
        struct_mat_id);

    return constituent;
  };

  if (scatra_material->material_type() == Core::Materials::m_matlist)
  {
    const auto actmat = std::dynamic_pointer_cast<const Mat::MatList>(scatra_material);
    for (int k = 0; k < my::numscal_; ++k)
      constituent_[k] = resolve(actmat->material_by_id(actmat->mat_id(k)));
  }
  else
  {
    constituent_[0] = resolve(scatra_material);
  }

  return 0;
}

template <Core::FE::CellType distype, int probdim>
double Discret::Elements::ScaTraEleCalcGrowthRemodel<distype,
    probdim>::get_stress_dependent_reaction_coeff(const Mixture::MixtureConstituentRemodelFiberSsi&
                                                      constituent,
    const Mat::PAR::ScatraGrowthRemodelMat::ScalarQuantity scalar_type, const int gp)
{
  if (gp == -1)
    FOUR_C_THROW(
        "EVALUATION_MAT: element_center is not implemented for ScaTraEleCalcGrowthRemodel");

  switch (scalar_type)
  {
    case Mat::PAR::ScatraGrowthRemodelMat::ScalarQuantity::growth:
      return -constituent.evaluate_growth_reaction_coefficient(gp);

    case Mat::PAR::ScatraGrowthRemodelMat::ScalarQuantity::remodeling:
      return -constituent.evaluate_remodeling_reaction_coefficient(gp);

    default:
      FOUR_C_THROW("Neither growth nor remodeling is activated for this material!");
  }
}

template <Core::FE::CellType distype, int probdim>
std::shared_ptr<Mat::Mixture>
Discret::Elements::ScaTraEleCalcGrowthRemodel<distype, probdim>::get_struct_material(
    const Core::Elements::Element* ele)
{
  if (ele->num_material() <= 1)
    FOUR_C_THROW(
        "ScaTra GR coupling requires a second structural material, but this is not defined for "
        "the element");

  auto mixture = std::dynamic_pointer_cast<Mat::Mixture>(ele->material(1));
  if (!mixture)
    FOUR_C_THROW(
        "ScaTra GR coupling expects material(1) to be a Mat::Mixture structural material used to "
        "query growth/remodeling reaction coefficients, but material(1) is not Mat::Mixture.");

  return mixture;
}


// template classes

// 1D elements
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::line2, 1>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::line2, 2>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::line2, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::tri3, 2>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::tri3, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::tri6, 2>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::quad4, 2>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::quad4, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::quad9, 2>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::hex8, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::hex27, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::tet4, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::tet10, 3>;
template class Discret::Elements::ScaTraEleCalcGrowthRemodel<Core::FE::CellType::pyramid5, 3>;

FOUR_C_NAMESPACE_CLOSE
