/*--------------------------------------------------------------------------*/
/*! \file

\brief Element evaluations for loma problems

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "4C_scatra_ele_calc_loma.hpp"

#include "4C_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_discretization_geometry_position_array.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_arrhenius_pv.hpp"
#include "4C_mat_arrhenius_spec.hpp"
#include "4C_mat_arrhenius_temp.hpp"
#include "4C_mat_ferech_pv.hpp"
#include "4C_mat_mixfrac.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_mat_tempdepwater.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_mat_yoghurt.hpp"
#include "4C_nurbs_discret_nurbs_utils.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_parameter_turbulence.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcLoma<distype>* DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcLoma<distype>>(
            new ScaTraEleCalcLoma<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::ScaTraEleCalcLoma(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
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
  if (my::turbparams_->MfsConservative())
    FOUR_C_THROW("Conservative formulation not supported for loma!");

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate single loma material  (protected)                vg 12/13  |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::materials(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  if (material->MaterialType() == CORE::Materials::m_mixfrac)
    mat_mix_frac(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_sutherland)
    mat_sutherland(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_tempdepwater)
    mat_temp_dep_water(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_arrhenius_pv)
    mat_arrhenius_pv(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_arrhenius_spec)
    mat_arrhenius_spec(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_arrhenius_temp)
    mat_arrhenius_temp(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_ferech_pv)
    mat_arrhenius_pv(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_thermostvenant)
    mat_thermo_st_venant_kirchhoff(material, k, densn, densnp, densam, visc);
  else if (material->MaterialType() == CORE::Materials::m_yoghurt)
    mat_yoghurt(material, k, densn, densnp, densam, visc);
  else
    FOUR_C_THROW("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 | material mixfrac                                            vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_mix_frac(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::MixFrac>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::MixFrac>(material);

  // compute mixture fraction at n+1 or n+alpha_F
  const double mixfracnp = my::scatravarmanager_->Phinp(0);

  // compute dynamic diffusivity at n+1 or n+alpha_F based on mixture fraction
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(mixfracnp), k);

  // compute density at n+1 or n+alpha_F based on mixture fraction
  densnp = actmat->ComputeDensity(mixfracnp);

  // set specific heat capacity at constant pressure to 1.0
  shc_ = 1.0;

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double mixfracam = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(mixfracam);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n
      const double mixfracn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(mixfracn);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp * densnp * actmat->EosFacA();

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(mixfracnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Sutherland                                         vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_sutherland(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::Sutherland>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = my::scatravarmanager_->Phinp(0);
  if (tempnp < 0.0) FOUR_C_THROW("Negative temperature in ScaTra Sutherland material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp), k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp, thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(tempam, thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(tempn, thermpressam_);
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
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material temperature-dependent water                        vg 07/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_temp_dep_water(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::TempDepWater>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::TempDepWater>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = my::scatravarmanager_->Phinp(0);
  if (tempnp < 0.0) FOUR_C_THROW("Negative temperature in ScaTra Sutherland material evaluation!");

  // compute diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp), k);

  // compute density at n+1 or n+alpha_F based on temperature
  densnp = actmat->ComputeDensity(tempnp);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(tempam);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(tempn);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Arrhenius PV                                       vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_arrhenius_pv(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::ArrheniusPV>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ArrheniusPV>(material);

  // get progress variable at n+1 or n+alpha_F
  const double provarnp = my::scatravarmanager_->Phinp(0);

  // get specific heat capacity at constant pressure and
  // compute temperature based on progress variable
  shc_ = actmat->ComputeShc(provarnp);

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = actmat->ComputeTemperature(provarnp);
  if (tempnp < 0.0)
    FOUR_C_THROW("Negative temperature in ScaTra Arrhenius progress-variable material evaluation!");

  // compute density at n+1 or n+alpha_F
  densnp = actmat->ComputeDensity(provarnp);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double provaram = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(provaram);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n
      const double provarn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(provarn);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp * actmat->ComputeFactor(provarnp);

  // compute diffusivity according to
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp), 0);

  // compute reaction coefficient for progress variable
  const double reacoef = actmat->compute_reaction_coeff(tempnp);

  // set different reaction terms in the reaction manager
  my::reamanager_->SetReaCoeff(reacoef, 0);

  // compute right-hand side contribution for progress variable
  // -> equal to reaction coefficient
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reacoef, 0);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Arrhenius Spec                                     vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_arrhenius_spec(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::ArrheniusSpec>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ArrheniusSpec>(material);

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = my::scatravarmanager_->Phinp(my::numscal_ - 1);
  if (tempnp < 0.0)
    FOUR_C_THROW("Negative temperature in ScaTra Arrhenius species material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp), k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp, thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[my::numscal_ - 1]);
    densam = actmat->ComputeDensity(tempam, thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(my::numscal_ - 1);
      densn = actmat->ComputeDensity(tempn, thermpressam_);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // factor for density gradient
  densgradfac_[k] = -densnp / tempnp;

  // compute reaction coefficient for species equation and set in reaction manager
  const double reacoef = actmat->compute_reaction_coeff(tempnp);
  my::reamanager_->SetReaCoeff(reacoef, k);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}

/*----------------------------------------------------------------------*
 | material Arrhenius Spec                                     vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_arrhenius_temp(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  if (k != my::numscal_ - 1)
    FOUR_C_THROW(
        "Temperature always needs to be the last variable for this Arrhenius-type system!");

  const Teuchos::RCP<const MAT::ArrheniusTemp>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ArrheniusTemp>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute species mass fraction and temperature at n+1 or n+alpha_F, including
  // check whether temperature is positive
  // (only two-equation systems, for the time being, such that only one species
  //  mass fraction possible)
  const double spmfnp = my::scatravarmanager_->Phinp(0);
  const double tempnp = my::scatravarmanager_->Phinp(k);
  if (tempnp < 0.0)
    FOUR_C_THROW("Negative temperature in ScaTra Arrhenius temperature material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp), k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp, thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[k]);
    densam = actmat->ComputeDensity(tempam, thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(k);
      densn = actmat->ComputeDensity(tempn, thermpressam_);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // factor for density gradient
  densgradfac_[k] = -densnp / tempnp;

  // compute sum of reaction rates for temperature equation divided by specific
  // heat capacity -> will be considered as a right-hand side contribution
  const double reatemprhs = actmat->ComputeReactionRHS(spmfnp, tempnp);
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reatemprhs, k);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Ferech PV                                          vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_ferech_pv(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::FerEchPV>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::FerEchPV>(material);

  // get progress variable at n+1 or n+alpha_F
  const double provarnp = my::scatravarmanager_->Phinp(0);

  // get specific heat capacity at constant pressure
  shc_ = actmat->ComputeShc(provarnp);

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = actmat->ComputeTemperature(provarnp);
  if (tempnp < 0.0)
    FOUR_C_THROW(
        "Negative temperature in ScaTra Ferziger and Echekki progress-variable material "
        "evaluation!");

  // compute density at n+1 or n+alpha_F
  densnp = actmat->ComputeDensity(provarnp);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double provaram = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(provaram);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n
      const double provarn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(provarn);
    }
    else
      densn = 1.0;
  }
  else
    densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp * actmat->ComputeFactor(provarnp);

  // compute diffusivity according to Ferech law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp), 0);

  // compute reaction coefficient for progress variable
  const double reacoef = actmat->compute_reaction_coeff(provarnp);

  // set different reaction terms in the reaction manager
  my::reamanager_->SetReaCoeff(reacoef, 0);

  // compute right-hand side contribution for progress variable
  // -> equal to reaction coefficient
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reacoef, 0);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material thermo St. Venant Kirchhoff                        vg 02/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_thermo_st_venant_kirchhoff(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::ThermoStVenantKirchhoff>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ThermoStVenantKirchhoff>(material);

  // get constant density
  densnp = actmat->Density();
  densam = densnp;
  densn = densnp;

  // set zero factor for density gradient
  densgradfac_[0] = 0.0;

  // set specific heat capacity at constant volume
  // (value divided by density here for its intended use on right-hand side)
  shc_ = actmat->Capacity() / densnp;

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
  my::reamanager_->SetReaCoeff(reacoef, 0);

  // ensure that temporal derivative of thermodynamic pressure is zero for
  // the present structure-based scalar transport
  thermpressdt_ = 0.0;

  // compute diffusivity as ratio of conductivity and specific heat capacity at constant volume
  my::diffmanager_->SetIsotropicDiff(actmat->Conductivity() / actmat->Capacity(), k);

  return;
}


/*----------------------------------------------------------------------*
 | material Yoghurt                                            vg 12/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::mat_yoghurt(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc     //!< fluid viscosity
)
{
  const Teuchos::RCP<const MAT::Yoghurt>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::Yoghurt>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(), 0);
  // diffus_[0] = actmat->ComputeDiffusivity();

  // get constant density
  densnp = actmat->Density();
  densam = densnp;
  densn = densnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or
      my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // compute temperature at n+1 or n+alpha_F and check whether it is positive
    const double tempnp = my::scatravarmanager_->Phinp(0);
    if (tempnp < 0.0) FOUR_C_THROW("Negative temperature in ScaTra yoghurt material evaluation!");

    // compute rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = this->get_strain_rate(my::evelnp_);

    // compute viscosity for Yoghurt-like flows according to Afonso et al. (2003)
    visc = actmat->ComputeViscosity(rateofstrain, tempnp);
  }

  return;
}


/*-----------------------------------------------------------------------------*
 | compute rhs containing bodyforce                                 ehrl 11/13 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::get_rhs_int(
    double& rhsint,       //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int k           //!< index of current scalar
)
{
  // get reatemprhs of species k from the reaction manager
  const double reatemprhs =
      Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->GetReaTempRhs(k);

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
      rhsint = my::bodyforce_[k].Dot(my::funct_) / shc_;
      rhsint += thermpressdt_ / shc_;
    }
    else
      rhsint = my::bodyforce_[k].Dot(my::funct_);
  }

  return;
}  // GetRhsInt


/*------------------------------------------------------------------------------------------*
 |  re-implementatio: calculation of convective element matrix: add conservative contributions  ehrl
 11/13 |
 *------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::calc_mat_conv_add_cons(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const CORE::LINALG::Matrix<nsd_, 1>& convelint, const CORE::LINALG::Matrix<nsd_, 1>& gradphi,
    const double vdiv, const double densnp, const double visc)
{
  // convective term using current scalar value
  const double cons_conv_phi = convelint.Dot(gradphi);

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
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::recompute_conv_phi_for_rhs(double& conv_phi,
    const int k, const CORE::LINALG::Matrix<nsd_, 1>& sgvelint,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const double densnp, const double densn,
    const double phinp, const double phin, const double vdiv)
{
  if (my::scatraparatimint_->IsIncremental())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->IsConservative())
    {
      // convective term in conservative form
      conv_phi += phinp * (vdiv + (densgradfac_[k] / densnp) * conv_phi);
    }

    // multiply convective term by density
    conv_phi *= densnp;
  }
  else if (not my::scatraparatimint_->IsIncremental() and my::scatraparatimint_->IsGenAlpha())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->IsConservative())
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
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
