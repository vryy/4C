/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation class containing routines for dissolving reactions within
       porous medium for ECM modeling

\level 2

 *----------------------------------------------------------------------*/

#include "4C_scatra_ele_calc_poro_reac_ECM.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_poro_ecm.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mat_structporo_reaction_ecm.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::ScaTraEleCalcPoroReacECM(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(
          numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(
          numdofpernode, numscal, disname)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>*
Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcPoroReacECM<distype>>(
            new ScaTraEleCalcPoroReacECM<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::materials(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  switch (material->MaterialType())
  {
    case Core::Materials::m_scatra:
      pororeac::mat_scatra(material, k, densn, densnp, densam, visc, iquad);
      break;
    default:
      FOUR_C_THROW("Material type %i is not supported", material->MaterialType());
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::get_material_params(
    const Core::Elements::Element* ele,  //!< the element we are dealing with
    std::vector<double>& densn,          //!< density at t_(n)
    std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,         //!< density at t_(n+alpha_M)
    double& visc,                        //!< fluid viscosity
    const int iquad                      //!< id of current gauss point
)
{
  // call poro base class to compute porosity
  poro::compute_porosity(ele);

  // get the material
  Teuchos::RCP<Core::Mat::Material> material = ele->Material();

  if (material->MaterialType() == Core::Materials::m_matlist_reactions)
  {
    const Teuchos::RCP<Mat::MatListReactions> actmat =
        Teuchos::rcp_dynamic_cast<Mat::MatListReactions>(material);
    if (actmat->NumMat() != my::numscal_) FOUR_C_THROW("Not enough materials in MatList.");

    for (int k = 0; k < actmat->NumReac(); ++k)
    {
      int matid = actmat->ReacID(k);
      Teuchos::RCP<Core::Mat::Material> singlemat = actmat->MaterialById(matid);

      Teuchos::RCP<Mat::ScatraMatPoroECM> scatramat =
          Teuchos::rcp_dynamic_cast<Mat::ScatraMatPoroECM>(singlemat);

      if (scatramat != Teuchos::null)
      {
        Teuchos::RCP<Mat::StructPoroReactionECM> structmat =
            Teuchos::rcp_dynamic_cast<Mat::StructPoroReactionECM>(my::ele_->Material(1));
        if (structmat == Teuchos::null) FOUR_C_THROW("cast to Mat::StructPoroReactionECM failed!");
        double structpot = compute_struct_chem_potential(structmat, iquad);

        scatramat->ComputeReacCoeff(structpot);
      }
    }
  }

  // call base class
  advreac::get_material_params(ele, densn, densnp, densam, visc, iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 19/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
double Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::compute_struct_chem_potential(
    Teuchos::RCP<Mat::StructPoroReactionECM>& structmat, const int gp)
{
  // gauss point displacements
  Core::LinAlg::Matrix<nsd_, 1> dispint(false);
  dispint.multiply(my::edispnp_, my::funct_);

  // transposed jacobian "dX/ds"
  Core::LinAlg::Matrix<nsd_, nsd_> xjm0;
  xjm0.multiply_nt(my::deriv_, poro::xyze0_);

  // inverse of transposed jacobian "ds/dX"
  Core::LinAlg::Matrix<nsd_, nsd_> xji0(true);
  xji0.invert(xjm0);

  // inverse of transposed jacobian "ds/dX"
  const double det0 = xjm0.determinant();

  my::xjm_.multiply_nt(my::deriv_, my::xyze_);
  const double det = my::xjm_.determinant();

  // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  const double J = det / det0;

  // ----------------------compute derivatives N_XYZ_ at gp w.r.t. material coordinates
  /// first derivatives of shape functions w.r.t. material coordinates
  Core::LinAlg::Matrix<nsd_, nen_> N_XYZ;
  N_XYZ.multiply(xji0, my::deriv_);

  // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ *
  // N_XYZ_^T
  static Core::LinAlg::Matrix<nsd_, nsd_> defgrd(false);
  defgrd.multiply_nt(my::xyze_, N_XYZ);

  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  static Core::LinAlg::Matrix<6, 1> glstrain(true);
  glstrain.clear();
  // if (kinemtype_ == Inpar::Solid::KinemType::nonlinearTotLag)
  {
    // Right Cauchy-Green tensor = F^T * F
    Core::LinAlg::Matrix<nsd_, nsd_> cauchygreen;
    cauchygreen.multiply_tn(defgrd, defgrd);
    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    if (nsd_ == 3)
    {
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);
    }
    else if (nsd_ == 2)
    {
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.0;
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = 0.0;
      glstrain(5) = 0.0;
    }
  }

  // fluid pressure at gauss point
  const double pres = my::eprenp_.dot(my::funct_);

  double pot = 0.0;

  structmat->ChemPotential(
      glstrain, poro::diff_manager()->GetPorosity(0), pres, J, my::eid_, pot, gp);

  return pot;
}


// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::quad4>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::quad9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::hex8>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::tet10>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::pyramid5>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::nurbs9>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcPoroReacECM<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
