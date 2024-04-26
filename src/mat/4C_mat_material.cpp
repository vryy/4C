/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_material.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_aaa_mixedeffects.hpp"
#include "4C_mat_aaagasser.hpp"
#include "4C_mat_aaaneohooke.hpp"
#include "4C_mat_aaaneohooke_stopro.hpp"
#include "4C_mat_aaaraghavanvorp_damage.hpp"
#include "4C_mat_air_0d_O2_saturation.hpp"
#include "4C_mat_arrhenius_pv.hpp"
#include "4C_mat_arrhenius_spec.hpp"
#include "4C_mat_arrhenius_temp.hpp"
#include "4C_mat_beam3r_plasticity.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"
#include "4C_mat_carreauyasuda.hpp"
#include "4C_mat_cnst_1d_art.hpp"
#include "4C_mat_constraintmixture.hpp"
#include "4C_mat_crosslinkermat.hpp"
#include "4C_mat_crystal_plasticity.hpp"
#include "4C_mat_damage.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_mat_electromagnetic.hpp"
#include "4C_mat_ferech_pv.hpp"
#include "4C_mat_fluid_linear_density_viscosity.hpp"
#include "4C_mat_fluid_murnaghantait.hpp"
#include "4C_mat_fluid_weakly_compressible.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_fluidporo_multiphase_reactions.hpp"
#include "4C_mat_fluidporo_multiphase_singlereaction.hpp"
#include "4C_mat_fluidporo_relpermeability_law.hpp"
#include "4C_mat_fluidporo_singlephase.hpp"
#include "4C_mat_fluidporo_singlephaseDof.hpp"
#include "4C_mat_fluidporo_singlephaselaw.hpp"
#include "4C_mat_fluidporo_viscosity_law.hpp"
#include "4C_mat_fourieriso.hpp"
#include "4C_mat_growth.hpp"
#include "4C_mat_growthremodel_elasthyper.hpp"
#include "4C_mat_hemoglobin_0d_O2_saturation.hpp"
#include "4C_mat_herschelbulkley.hpp"
#include "4C_mat_ion.hpp"
#include "4C_mat_lin_elast_1D.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_list_chemoreac.hpp"
#include "4C_mat_list_chemotaxis.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_lubrication_law.hpp"
#include "4C_mat_lubrication_mat.hpp"
#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_mat_maxwell_0d_acinus_DoubleExponential.hpp"
#include "4C_mat_maxwell_0d_acinus_Exponential.hpp"
#include "4C_mat_maxwell_0d_acinus_NeoHookean.hpp"
#include "4C_mat_maxwell_0d_acinus_Ogden.hpp"
#include "4C_mat_membrane_active_strain.hpp"
#include "4C_mat_membrane_elasthyper.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_mat_mixfrac.hpp"
#include "4C_mat_mixture.hpp"
#include "4C_mat_modpowerlaw.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_muscle_combo.hpp"
#include "4C_mat_muscle_giantesio.hpp"
#include "4C_mat_muscle_weickenmeier.hpp"
#include "4C_mat_myocard.hpp"
#include "4C_mat_newman.hpp"
#include "4C_mat_newman_multiscale.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_par_parameter.hpp"
#include "4C_mat_particle_dem.hpp"
#include "4C_mat_particle_sph_boundary.hpp"
#include "4C_mat_particle_sph_fluid.hpp"
#include "4C_mat_particle_wall_dem.hpp"
#include "4C_mat_permeablefluid.hpp"
#include "4C_mat_plastic_VarConstUpdate.hpp"
#include "4C_mat_plasticdruckerprager.hpp"
#include "4C_mat_plasticelasthyper.hpp"
#include "4C_mat_plasticlinelast.hpp"
#include "4C_mat_plasticnlnlogneohooke.hpp"
#include "4C_mat_poro_density_law.hpp"
#include "4C_mat_poro_law.hpp"
#include "4C_mat_robinson.hpp"
#include "4C_mat_scalardepinterp.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_aniso.hpp"
#include "4C_mat_scatra_chemotaxis.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_mat_scatra_multiscale.hpp"
#include "4C_mat_scatra_poro_ecm.hpp"
#include "4C_mat_scatra_reaction.hpp"
#include "4C_mat_scl.hpp"
#include "4C_mat_soret.hpp"
#include "4C_mat_spring.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mat_structporo_reaction.hpp"
#include "4C_mat_structporo_reaction_ecm.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_mat_superelastic_sma.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_mat_tempdepwater.hpp"
#include "4C_mat_thermoplastichyperelast.hpp"
#include "4C_mat_thermoplasticlinelast.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_mat_viscoanisotropic.hpp"
#include "4C_mat_viscoelasthyper.hpp"
#include "4C_mat_visconeohooke.hpp"
#include "4C_mat_viscoplastic_no_yield_surface.hpp"
#include "4C_mat_yoghurt.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::Material::Factory(int matnum)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case CORE::Materials::m_fluid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::NewtonianFluid(curmat));
      auto* params = static_cast<MAT::PAR::NewtonianFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluid_murnaghantait:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MurnaghanTaitFluid(curmat));
      auto* params = static_cast<MAT::PAR::MurnaghanTaitFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluid_linear_density_viscosity:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LinearDensityViscosity(curmat));
      auto* params = static_cast<MAT::PAR::LinearDensityViscosity*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluid_weakly_compressible:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::WeaklyCompressibleFluid(curmat));
      auto* params = static_cast<MAT::PAR::WeaklyCompressibleFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_stvenant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::StVenantKirchhoff(curmat));
      auto* params = static_cast<MAT::PAR::StVenantKirchhoff*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_thermostvenant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoStVenantKirchhoff(curmat));
      auto* params = static_cast<MAT::PAR::ThermoStVenantKirchhoff*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_thermopllinelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoPlasticLinElast(curmat));
      auto* params = static_cast<MAT::PAR::ThermoPlasticLinElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_pldruckprag:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticDruckerPrager(curmat));
      auto* params = static_cast<MAT::PAR::PlasticDruckerPrager*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_thermoplhyperelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoPlasticHyperElast(curmat));
      auto* params = static_cast<MAT::PAR::ThermoPlasticHyperElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_plnlnlogneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticNlnLogNeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::PlasticNlnLogNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_pllinelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticLinElast(curmat));
      auto* params = static_cast<MAT::PAR::PlasticLinElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_vp_no_yield_surface:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ViscoPlasticNoYieldSurface(curmat));
      auto* params = static_cast<MAT::PAR::ViscoPlasticNoYieldSurface*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_vp_robinson:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Robinson(curmat));
      auto* params = static_cast<MAT::PAR::Robinson*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_elpldamage:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Damage(curmat));
      auto* params = static_cast<MAT::PAR::Damage*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_struct_multiscale:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MicroMaterial(curmat));
      auto* params = static_cast<MAT::PAR::MicroMaterial*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_visconeohooke:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ViscoNeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::ViscoNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_viscoanisotropic:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ViscoAnisotropic(curmat));
      auto* params = static_cast<MAT::PAR::ViscoAnisotropic*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_aaaneohooke:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::AAAneohooke(curmat));
      auto* params = static_cast<MAT::PAR::AAAneohooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_aaaneohooke_stopro:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::AaAneohookeStopro(curmat));
      auto* params = static_cast<MAT::PAR::AaAneohookeStopro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_aaagasser:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::AAAgasser(curmat));
      auto* params = static_cast<MAT::PAR::AAAgasser*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_aaaraghavanvorp_damage:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::AaAraghavanvorpDamage(curmat));
      auto* params = static_cast<MAT::PAR::AaAraghavanvorpDamage*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_aaa_mixedeffects:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::AaaMixedeffects(curmat));
      auto* params = static_cast<MAT::PAR::AaaMixedeffects*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_lubrication:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationMat(curmat));
      auto* params = static_cast<MAT::PAR::LubricationMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_lubrication_law_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::LubricationLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_lubrication_law_barus:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawBarus(curmat));
      auto* params = static_cast<MAT::PAR::LubricationLawBarus*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_lubrication_law_roeland:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawRoeland(curmat));
      auto* params = static_cast<MAT::PAR::LubricationLawRoeland*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ScatraMat(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_reaction_poroECM:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatPoroECM(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatPoroECM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_reaction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraReactionMat(curmat));
      auto* params = static_cast<MAT::PAR::ScatraReactionMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_multiporo_fluid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroFluid(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_multiporo_volfrac:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroVolFrac(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroVolFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_multiporo_solid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroSolid(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroSolid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_multiporo_temperature:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroTemperature(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroTemperature*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_multiscale:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMultiScale(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMultiScale*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_chemotaxis:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraChemotaxisMat(curmat));
      auto* params = static_cast<MAT::PAR::ScatraChemotaxisMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scatra_aniso:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatAniso(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatAniso*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_muscle_combo:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MuscleCombo(curmat));
      auto* params = static_cast<MAT::PAR::MuscleCombo*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_muscle_giantesio:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MuscleGiantesio(curmat));
      auto* params = static_cast<MAT::PAR::MuscleGiantesio*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_muscle_weickenmeier:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MuscleWeickenmeier(curmat));
      auto* params = static_cast<MAT::PAR::MuscleWeickenmeier*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_myocard:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Myocard(curmat));
      auto* params = static_cast<MAT::PAR::Myocard*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_mixfrac:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MixFrac(curmat));
      auto* params = static_cast<MAT::PAR::MixFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_sutherland:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Sutherland(curmat));
      auto* params = static_cast<MAT::PAR::Sutherland*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_tempdepwater:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::TempDepWater(curmat));
      auto* params = static_cast<MAT::PAR::TempDepWater*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_arrhenius_spec:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ArrheniusSpec(curmat));
      auto* params = static_cast<MAT::PAR::ArrheniusSpec*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_arrhenius_temp:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ArrheniusTemp(curmat));
      auto* params = static_cast<MAT::PAR::ArrheniusTemp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_arrhenius_pv:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ArrheniusPV(curmat));
      auto* params = static_cast<MAT::PAR::ArrheniusPV*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_ferech_pv:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FerEchPV(curmat));
      auto* params = static_cast<MAT::PAR::FerEchPV*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_carreauyasuda:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::CarreauYasuda(curmat));
      auto* params = static_cast<MAT::PAR::CarreauYasuda*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_modpowerlaw:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ModPowerLaw(curmat));
      auto* params = static_cast<MAT::PAR::ModPowerLaw*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_herschelbulkley:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::HerschelBulkley(curmat));
      auto* params = static_cast<MAT::PAR::HerschelBulkley*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_yoghurt:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Yoghurt(curmat));
      auto* params = static_cast<MAT::PAR::Yoghurt*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_permeable_fluid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PermeableFluid(curmat));
      auto* params = static_cast<MAT::PAR::PermeableFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FluidPoro(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_multiphase:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroMultiPhase(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroMultiPhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_multiphase_reactions:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroMultiPhaseReactions(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroMultiPhaseReactions*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_singlereaction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroSingleReaction(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroSingleReaction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_singlephase:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroSinglePhase(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroSinglePhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_singlevolfrac:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroSingleVolFrac(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroSingleVolFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_volfracpressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroVolFracPressure(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroVolFracPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_law_linear:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::PoroLawLinear(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawLinear*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_law_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_law_logNeoHooke_Penalty:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawNeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_law_incompr_skeleton:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawIncompSkeleton(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawIncompSkeleton*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_law_linear_biot:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawLinBiot(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawLinBiot*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_law_density_dependent:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawDensityDependent(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawDensityDependent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_densitylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::PoroDensityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_poro_densitylaw_exp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawExp(curmat));
      auto* params = static_cast<MAT::PAR::PoroDensityLawExp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phaselaw_linear:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawLinear(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawLinear*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phaselaw_tangent:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawTangent(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawTangent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phaselaw_constraint:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawConstraint(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawConstraint*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phaselaw_byfunction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawByFunction(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawByFunction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_relpermeabilitylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawConstant(curmat));
      auto* params =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_relpermeabilitylaw_exp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawExponent(curmat));
      auto* params =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawExponent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_viscositylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroViscosityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_viscositylaw_celladh:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawCellAdherence(curmat));
      auto* params =
          static_cast<MAT::PAR::FluidPoroViscosityLawCellAdherence*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phasedof_diffpressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofDiffPressure(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseDofDiffPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phasedof_pressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressure(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseDofPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_fluidporo_phasedof_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofSaturation(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseDofSaturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_matlist:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MatList(curmat));
      auto* params = static_cast<MAT::PAR::MatList*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_matlist_reactions:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MatListReactions(curmat));
      auto* params = dynamic_cast<MAT::PAR::MatListReactions*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_matlist_chemotaxis:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MatListChemotaxis(curmat));
      auto* params = dynamic_cast<MAT::PAR::MatListChemotaxis*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_matlist_chemoreac:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MatListChemoReac(curmat));
      auto* params = dynamic_cast<MAT::PAR::MatListChemoReac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_elchmat:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ElchMat(curmat));
      auto* params = static_cast<MAT::PAR::ElchMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_elchphase:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ElchPhase(curmat));
      auto* params = static_cast<MAT::PAR::ElchPhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_ion:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Ion(curmat));
      auto* params = static_cast<MAT::PAR::Ion*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_electrode:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Electrode(curmat));
      auto* params = static_cast<MAT::PAR::Electrode*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_newman:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Newman(curmat));
      auto* params = static_cast<MAT::PAR::Newman*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_newman_multiscale:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::NewmanMultiScale(curmat));
      auto* params = static_cast<MAT::PAR::NewmanMultiScale*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_scl:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Scl(curmat));
      auto* params = static_cast<MAT::PAR::Scl*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_elasthyper:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_viscoelasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ViscoElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::ViscoElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_plelasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::PlasticElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_plelasthyperVCU:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticElastHyperVCU(curmat));
      auto* params = static_cast<MAT::PAR::PlasticElastHyperVCU*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_sc_dep_interp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScalarDepInterp(curmat));
      auto* params = static_cast<MAT::PAR::ScalarDepInterp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_structporo:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::StructPoro(curmat));
      auto* params = static_cast<MAT::PAR::StructPoro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_structpororeaction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::StructPoroReaction(curmat));
      auto* params = static_cast<MAT::PAR::StructPoroReaction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_structpororeactionECM:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::StructPoroReactionECM(curmat));
      auto* params = static_cast<MAT::PAR::StructPoroReactionECM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::mes_couplogneohooke:
    case CORE::Materials::mes_couplogmixneohooke:
    case CORE::Materials::mes_coupexppol:
    case CORE::Materials::mes_coupneohooke:
    case CORE::Materials::mes_coupblatzko:
    case CORE::Materials::mes_isoneohooke:
    case CORE::Materials::mes_isoogden:
    case CORE::Materials::mes_isoyeoh:
    case CORE::Materials::mes_iso1pow:
    case CORE::Materials::mes_iso2pow:
    case CORE::Materials::mes_coup1pow:
    case CORE::Materials::mes_coup2pow:
    case CORE::Materials::mes_coup3pow:
    case CORE::Materials::mes_coup13apow:
    case CORE::Materials::mes_coupmooneyrivlin:
    case CORE::Materials::mes_coupsimopister:
    case CORE::Materials::mes_isoexpopow:
    case CORE::Materials::mes_isomooneyrivlin:
    case CORE::Materials::mes_isomuscleblemker:
    case CORE::Materials::mes_volsussmanbathe:
    case CORE::Materials::mes_volpenalty:
    case CORE::Materials::mes_vologden:
    case CORE::Materials::mes_volpow:
    case CORE::Materials::mes_anisoactivestress_evolution:
    case CORE::Materials::mes_coupanisoexpoactive:
    case CORE::Materials::mes_coupanisoexpo:
    case CORE::Materials::mes_coupanisoexposhear:
    case CORE::Materials::mes_coupanisopow:
    case CORE::Materials::mes_couptransverselyisotropic:
    case CORE::Materials::mes_coupanisoexpotwocoup:
    case CORE::Materials::mes_coupanisoneohooke:
    case CORE::Materials::mes_coupanisoneohooke_varprop:
    case CORE::Materials::mes_isoanisoexpo:
    case CORE::Materials::mes_structuraltensorstratgy:
    case CORE::Materials::mes_coupvarga:
    case CORE::Materials::mes_isovarga:
    case CORE::Materials::mes_isovolaaagasser:
    case CORE::Materials::mes_isotestmaterial:
    case CORE::Materials::mes_coupmyocard:
    case CORE::Materials::mes_isoratedep:
    case CORE::Materials::mes_genmax:
    case CORE::Materials::mes_fract:
    case CORE::Materials::mes_generalizedgenmax:
    case CORE::Materials::mes_viscopart:
    case CORE::Materials::mes_viscobranch:
    case CORE::Materials::mes_remodelfiber:
    case CORE::Materials::m_growth_aniso_strain:
    case CORE::Materials::m_growth_aniso_stress:
    case CORE::Materials::m_growth_aniso_strain_const_trig:
    case CORE::Materials::m_growth_aniso_stress_const_trig:
    case CORE::Materials::m_growth_iso_stress:
    case CORE::Materials::m_growth_ac:
    case CORE::Materials::m_growth_ac_radial:
    case CORE::Materials::m_growth_ac_radial_refconc:
    case CORE::Materials::m_growth_const:
    case CORE::Materials::mes_coupSVK:
    case CORE::Materials::mfi_lin_scalar_aniso:
    case CORE::Materials::mfi_lin_scalar_iso:
    case CORE::Materials::mfi_lin_temp_iso:
    case CORE::Materials::mfi_no_growth:
    case CORE::Materials::mfi_time_funct:
    case CORE::Materials::mfi_poly_intercal_frac_aniso:
    case CORE::Materials::mfi_poly_intercal_frac_iso:
    case CORE::Materials::mix_rule_function:
    case CORE::Materials::mix_rule_map:
    case CORE::Materials::mix_rule_simple:
    case CORE::Materials::mix_rule_growthremodel:
    case CORE::Materials::mix_elasthyper:
    case CORE::Materials::mix_elasthyper_damage:
    case CORE::Materials::mix_elasthyper_elastin_membrane:
    case CORE::Materials::mix_full_constrained_mixture_fiber:
    case CORE::Materials::mix_solid_material:
    case CORE::Materials::mix_growth_strategy_anisotropic:
    case CORE::Materials::mix_growth_strategy_isotropic:
    case CORE::Materials::mix_growth_strategy_stiffness:
    case CORE::Materials::mix_prestress_strategy_constant:
    case CORE::Materials::mix_prestress_strategy_cylinder:
    case CORE::Materials::mix_prestress_strategy_iterative:
    case CORE::Materials::mix_remodelfiber_expl:
    case CORE::Materials::mix_remodelfiber_impl:
    case CORE::Materials::mix_remodelfiber_material_exponential:
    case CORE::Materials::mix_remodelfiber_material_exponential_active:
    {
      return Teuchos::null;
    }
    case CORE::Materials::m_cnst_art:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Cnst1dArt(curmat));
      auto* params = static_cast<MAT::PAR::Cnst1dArt*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_maxwell_acinus:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell0dAcinus(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell0dAcinus*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_maxwell_acinus_neohookean:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell0dAcinusNeoHookean(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell0dAcinusNeoHookean*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_maxwell_acinus_exponential:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell0dAcinusExponential(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell0dAcinusExponential*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_maxwell_acinus_doubleexponential:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell0dAcinusDoubleExponential(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell0dAcinusDoubleExponential*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_maxwell_acinus_ogden:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell0dAcinusOgden(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell0dAcinusOgden*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_o2_hemoglobin_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Hemoglobin0dO2Saturation(curmat));
      auto* params = static_cast<MAT::PAR::Hemoglobin0dO2Saturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_0d_o2_air_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Air0dO2Saturation(curmat));
      auto* params = static_cast<MAT::PAR::Air0dO2Saturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_th_fourier_iso:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FourierIso(curmat));
      auto* params = static_cast<MAT::PAR::FourierIso*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_soret:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Soret(curmat));
      auto* params = static_cast<MAT::PAR::Soret*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_membrane_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MembraneElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::MembraneElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_membrane_activestrain:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MembraneActiveStrain(curmat));
      auto* params = static_cast<MAT::PAR::MembraneActiveStrain*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_growthremodel_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthRemodelElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::GrowthRemodelElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_mixture:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Mixture(curmat));
      auto* params = dynamic_cast<MAT::PAR::Mixture*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_multiplicative_split_defgrad_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MultiplicativeSplitDefgradElastHyper(curmat));
      auto* params =
          static_cast<MAT::PAR::MultiplicativeSplitDefgradElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_growth_volumetric:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Growth(curmat));
      auto* params = static_cast<MAT::PAR::Growth*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_constraintmixture:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ConstraintMixture(curmat));
      auto* params = static_cast<MAT::PAR::ConstraintMixture*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_beam_reissner_elast_hyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastHyperMaterialParams(curmat));
      auto* params =
          static_cast<MAT::PAR::BeamReissnerElastHyperMaterialParams*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_beam_reissner_elast_plastic:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastPlasticMaterialParams(curmat));
      auto* params =
          static_cast<MAT::PAR::BeamReissnerElastPlasticMaterialParams*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_beam_reissner_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode(curmat));
      auto* params_bymode =
          static_cast<MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode*>(curmat->Parameter());
      return params_bymode->CreateMaterial();
    }
    case CORE::Materials::m_beam_kirchhoff_elast_hyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamKirchhoffElastHyperMaterialParams(curmat));
      auto* params =
          static_cast<MAT::PAR::BeamKirchhoffElastHyperMaterialParams*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_beam_kirchhoff_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode(curmat));
      auto* params_bymode =
          static_cast<MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode*>(curmat->Parameter());
      return params_bymode->CreateMaterial();
    }
    case CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(
            new MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams(curmat));
      auto* params = static_cast<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams*>(
          curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(
            new MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(curmat));
      auto* params_bymode =
          static_cast<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode*>(
              curmat->Parameter());
      return params_bymode->CreateMaterial();
    }
    case CORE::Materials::m_crosslinkermat:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::CrosslinkerMat(curmat));
      auto* params = static_cast<MAT::PAR::CrosslinkerMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_spring:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Spring(curmat));
      auto* params = static_cast<MAT::PAR::Spring*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_particle_sph_fluid:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialSPHFluid(curmat));
      auto* params = dynamic_cast<MAT::PAR::ParticleMaterialSPHFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_particle_sph_boundary:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialSPHBoundary(curmat));
      auto* params = dynamic_cast<MAT::PAR::ParticleMaterialSPHBoundary*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_particle_dem:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialDEM(curmat));
      auto* params = dynamic_cast<MAT::PAR::ParticleMaterialDEM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_particle_wall_dem:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleWallMaterialDEM(curmat));
      auto* params = static_cast<MAT::PAR::ParticleWallMaterialDEM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_electromagneticmat:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ElectromagneticMat(curmat));
      auto* params = static_cast<MAT::PAR::ElectromagneticMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_superelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::SuperElasticSMA(curmat));
      auto* params = static_cast<MAT::PAR::SuperElasticSMA*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_crystplast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::CrystalPlasticity(curmat));
      auto* params = dynamic_cast<MAT::PAR::CrystalPlasticity*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_linelast1D:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::LinElast1D(curmat));
      auto* params = static_cast<MAT::PAR::LinElast1D*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case CORE::Materials::m_linelast1D_growth:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LinElast1DGrowth(curmat));
      auto* params = static_cast<MAT::PAR::LinElast1DGrowth*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    default:
      FOUR_C_THROW("unknown material type %d", curmat->Type());
      break;
  }

  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
