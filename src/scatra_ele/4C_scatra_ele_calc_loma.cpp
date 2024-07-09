/*--------------------------------------------------------------------------*/
/*! \file

\brief Element evaluations for loma problems

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "4C_scatra_ele_calc_loma.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_fem_nurbs_discretization_utils.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_parameter_turbulence.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcLoma<distype>*
Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcLoma<distype>>(
            new ScaTraEleCalcLoma<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::ScaTraEleCalcLoma(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      ephiam_(my::numscal_),
      densgradfac_(my::numscal_, 0.0),
      thermpressnp_(0.0),
      thermpressam_(0.0),
      thermpressdt_(0.0),
      shc_(1.0)
{
  // set appropriate reaction manager
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerLoma(my::numscal_));

  // safety check
  if (my::turbparams_->mfs_conservative())
    FOUR_C_THROW("Conservative formulation not supported for loma!");

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate single loma material  (protected)                vg 12/13  |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::materials(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  if (material->material_type() == Core::Materials::m_sutherland)
    mat_sutherland(material, k, densn, densnp, densam, visc);
  else if (material->material_type() == Core::Materials::m_thermostvenant)
    mat_thermo_st_venant_kirchhoff(material, k, densn, densnp, densam, visc);
  else
    FOUR_C_THROW("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 | material Sutherland                                         vg 12/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_sutherland(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const Mat::Sutherland>& actmat =
      Teuchos::rcp_dynamic_cast<const Mat::Sutherland>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->shc();

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = my::scatravarmanager_->phinp(0);
  if (tempnp < 0.0) FOUR_C_THROW("Negative temperature in ScaTra Sutherland material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->set_isotropic_diff(actmat->compute_diffusivity(tempnp), k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->compute_density(tempnp, thermpressnp_);

  if (my::scatraparatimint_->is_gen_alpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.dot(ephiam_[0]);
    densam = actmat->compute_density(tempam, thermpressam_);

    if (not my::scatraparatimint_->is_incremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->phin(0);
      densn = actmat->compute_density(tempn, thermpressam_);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp / tempnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->rb_sub_gr_vel() or
      my::turbparams_->turb_model() == Inpar::FLUID::multifractal_subgrid_scales)
    visc = actmat->compute_viscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material thermo St. Venant Kirchhoff                        vg 02/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_thermo_st_venant_kirchhoff(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const Mat::ThermoStVenantKirchhoff>& actmat =
      Teuchos::rcp_dynamic_cast<const Mat::ThermoStVenantKirchhoff>(material);

  // get constant density
  densnp = actmat->density();
  densam = densnp;
  densn = densnp;

  // set zero factor for density gradient
  densgradfac_[0] = 0.0;

  // set specific heat capacity at constant volume
  // (value divided by density here for its intended use on right-hand side)
  shc_ = actmat->capacity() / densnp;

  // compute velocity divergence required for reaction coefficient
  // double vdiv(0.0);
  // get_divergence(vdiv,evelnp_);

  // compute reaction coefficient
  // (divided by density due to later multiplication by density in CalMatAndRHS)
  // const double reacoef = -vdiv_*actmat->st_modulus()/(actmat->Capacity()*densnp);
  const double reacoef = 0.0;

  // set reaction flag to true, check whether reaction coefficient is positive
  // and set derivative of reaction coefficient
  // if (reacoef > 1e-14) reaction_ = true;
  // if (reacoef < -1e-14)
  //  FOUR_C_THROW("Reaction coefficient for Thermo St. Venant-Kirchhoff material is not positive:
  //  %f",0, reacoef);
  // reacoeffderiv_[0] = reacoef;

  // set different reaction terms in the reaction manager
  my::reamanager_->set_rea_coeff(reacoef, 0);

  // ensure that temporal derivative of thermodynamic pressure is zero for
  // the present structure-based scalar transport
  thermpressdt_ = 0.0;

  // compute diffusivity as ratio of conductivity and specific heat capacity at constant volume
  my::diffmanager_->set_isotropic_diff(actmat->conductivity() / actmat->capacity(), k);

  return;
}


/*-----------------------------------------------------------------------------*
 | compute rhs containing bodyforce                                 ehrl 11/13 |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::get_rhs_int(
    double& rhsint,       //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int k           //!< index of current scalar
)
{
  // get reatemprhs of species k from the reaction manager
  const double reatemprhs =
      Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->get_rea_temp_rhs(k);

  // Three cases have to be distinguished for computing the rhs:
  // 1) reactive temperature equation: reaction-rate term
  //    (divided by specific heat capacity)
  // 2) non-reactive temperature equation: heat-source term and
  //    temporal derivative of thermodynamic pressure
  //    (both divided by specific heat capacity)
  // 3) species equation: only potential body force (usually zero)
  const double tol = 1e-8;
  if ((reatemprhs < (0.0 - tol)) or (reatemprhs > (0.0 + tol)))
    rhsint = densnp * reatemprhs / shc_;
  else
  {
    if (k == my::numscal_ - 1)
    {
      rhsint = my::bodyforce_[k].dot(my::funct_) / shc_;
      rhsint += thermpressdt_ / shc_;
    }
    else
      rhsint = my::bodyforce_[k].dot(my::funct_);
  }

  return;
}  // GetRhsInt


/*------------------------------------------------------------------------------------------*
 |  re-implementatio: calculation of convective element matrix: add conservative contributions  ehrl
 11/13 |
 *------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::calc_mat_conv_add_cons(
    Core::LinAlg::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const Core::LinAlg::Matrix<nsd_, 1>& convelint, const Core::LinAlg::Matrix<nsd_, 1>& gradphi,
    const double vdiv, const double densnp, const double visc)
{
  // convective term using current scalar value
  const double cons_conv_phi = convelint.dot(gradphi);

  const double consfac = timefacfac * (densnp * vdiv + densgradfac_[k] * cons_conv_phi);
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = consfac * my::funct_(vi);
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + k;

      emat(fvi, fui) += v * my::funct_(ui);
    }
  }
  return;
}


/*------------------------------------------------------------------- *
 | re-implementatio: adaption of convective term for rhs   ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::recompute_conv_phi_for_rhs(double& conv_phi,
    const int k, const Core::LinAlg::Matrix<nsd_, 1>& sgvelint,
    const Core::LinAlg::Matrix<nsd_, 1>& gradphi, const double densnp, const double densn,
    const double phinp, const double phin, const double vdiv)
{
  if (my::scatraparatimint_->is_incremental())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->is_conservative())
    {
      // convective term in conservative form
      conv_phi += phinp * (vdiv + (densgradfac_[k] / densnp) * conv_phi);
    }

    // multiply convective term by density
    conv_phi *= densnp;
  }
  else if (not my::scatraparatimint_->is_incremental() and my::scatraparatimint_->is_gen_alpha())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->is_conservative())
    {
      // convective term in conservative form
      // caution: velocity divergence is for n+1 and not for n!
      // -> hopefully, this inconsistency is of small amount
      conv_phi += phin * (vdiv + (densgradfac_[k] / densnp) * conv_phi);
    }

    // multiply convective term by density
    conv_phi *= densn;
  }

  return;
}


// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::quad4>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::hex8>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tet10>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::pyramid5>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
